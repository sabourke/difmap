/* FITS Table extension HDU code */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "sysfits.h"
#include "fits.h"
#include "utils.h"
#include "ahdu.h"
#include "fitkey.h"

static int get_format(char *tform, Afield *field);
static long acol_return(char *tmpbuf, long iret);
static Afield *new_Afields(Ahdu *ahdu);

static COL_FINDFN(acol_find);
static COL_TYPEFN(acol_type);
static COL_DIMFN(acol_dim);
static COL_VALFN(acol_value);
static COL_NAMEFN(acol_name);
static COL_SETFN(acol_set);

Tabfn atabfn={acol_value, acol_find, acol_type, acol_dim, acol_name, acol_set};

static GETFN(get_ahdu);
static NEWFN(new_ahdu);
static DELFN(del_ahdu);
static SIZEFN(siz_ahdu);
static ADDFN(add_ahdu);
static COPYFN(cop_ahdu);
static ENDFN(end_ahdu);

Hdutab ahdufns={"TABLE", "TABLE", get_ahdu, new_ahdu, del_ahdu, siz_ahdu,
		add_ahdu, cop_ahdu, end_ahdu};

/* Create a scratch buffer to read and write column items via */

enum {MAX_ITEM_WIDTH=132};
static char item_buff[MAX_ITEM_WIDTH+1];

/*.......................................................................
 * Clear the derived part of an ASCII-table header descriptor to the
 * point at which it can be safely sent to del_Hdu().
 *
 * Input:
 *  hdu     Hdu *  The base-class pointer to a Phdu descriptor.
 * Output:
 *  return  Hdu *  'hdu', or NULL if an error occurs and hdu is deleted.
 */
static NEWFN(new_ahdu)
{
  Ahdu *ahdu = (Ahdu *) hdu;
  (void) new_table(hdu);
  ahdu->fields = NULL;
/*
 * Unlike other HDUs, ASCII tables must be padded with blanks (ASCII 32).
 */
  hdu->pad = 32;
  return hdu;
}

/*.......................................................................
 * Delete the derived parts of an ASCII-table HDU descriptor.
 *
 * Input:
 *  hdu      Hdu *  The base-class pointer to the Ahdu descriptor.
 */
static DELFN(del_ahdu)
{
  Ahdu *ahdu = (Ahdu *) hdu;
  Afield *field;
  int i;
/*
 * Delete the generic table members.
 */
  del_table(hdu);
/*
 * Clean-up the array of table field descriptors.
 */
  if(ahdu->fields) {
    for(i=0; i<ahdu->tfields; i++) {
      field = &ahdu->fields[i];
      if(field->tnull)
	free(field->tnull);
      if(field->tform)
	free(field->tform);
      if(field->ttype)
	free(field->ttype);
      if(field->tunit)
	free(field->tunit);
    };
    free(ahdu->fields);
  };
  return;
}

/*.......................................................................
 * Read the ASCII-table parts of a FITS header and record
 * the table description in an ASCII-table-HDU descriptor.
 *
 * Input:
 *  fits  Fits *  The descriptor of the FITS file.
 * Input/Output:
 *  hdu    Hdu *  The base-class descriptor to a Ahdu descriptor.
 * Output:
 *  return int    0 - OK.
 */
static GETFN(get_ahdu)
{
  Ahdu *ahdu = (Ahdu *) hdu;
  Afield *field;     /* A field descriptor */
  Fitkey key;        /* Keyword-value pair temporary storage */
  int i;
  /* Table field descriptor keywords */
  enum {AHDU_TBCOL, AHDU_TFORM, AHDU_TSCAL, AHDU_TZERO, AHDU_TNULL,
	AHDU_TTYPE, AHDU_TUNIT};
  Fitkey tkeys[]={
    {"TBCOL", 0, AHDU_TBCOL, DAT_INT, NULL, NULL},
    {"TFORM", 0, AHDU_TFORM, DAT_STR, NULL, NULL},
    {"TSCAL", 0, AHDU_TSCAL, DAT_DBL, NULL, NULL},
    {"TZERO", 0, AHDU_TZERO, DAT_DBL, NULL, NULL},
    {"TNULL", 0, AHDU_TNULL, DAT_STR, NULL, NULL},
    {"TTYPE", 0, AHDU_TTYPE, DAT_STR, NULL, NULL},
    {"TUNIT", 0, AHDU_TUNIT, DAT_STR, NULL, NULL},
  };
/*
 * Check the input Hdu type.
 */
  if(hdu->type != F_TABLE) {
    fprintf(stderr, "get_ahdu: Incompatible HDU descriptor received\n");
    return 1;
  };
/*
 * Check that NAXIS has the expected value.
 */
  if(hdu->naxis != 2) {
    fprintf(stderr, "Invalid NAXIS value (%d) != 2 in an ASCII table\n",
	    hdu->naxis);
    return 1;
  };
/*
 * The next two keywords should be PCOUNT and GCOUNT. These have
 * already been read by get_Hdu() so simply warn if they are not there.
 */
  if(get_key(fits, hdu, "PCOUNT", DAT_INT, NO_SEEK, &key))
    fprintf(stderr, "Missing PCOUNT keyword in FITS table header\n");
  if(get_key(fits, hdu, "GCOUNT", DAT_INT, NO_SEEK, &key))
    fprintf(stderr, "Missing GCOUNT keyword in FITS table header\n");
/*
 * Check the sanity of the PCOUNT and GCOUNT values.
 */
  if(ahdu->gcount != 1 || ahdu->pcount != 0) {
    fprintf(stderr, "Illegal values of PCOUNT=%d GCOUNT=%d in ascii table\n",
	    ahdu->pcount, ahdu->gcount);
    return 1;
  };
/*
 * Get the TFIELDS keyword (number of table fields per row).
 */
  if(get_key(fits, hdu, "TFIELDS", DAT_INT, NO_SEEK, &key)) {
    fprintf(stderr, "Missing TFIELDS keyword in FITS table header\n");
  } else {
    ahdu->tfields = KEYINT(key);
  };
/*
 * Allocate the array of field descriptors.
 */
  if(new_Afields(ahdu)==NULL)
    return 1;
/*
 * Get the table-field descriptor keywords. The TBCOL and TFORM keywords
 * are the only mandatory keywords.
 */
  if(ahdu->tfields > 0) {
    while(next_key(fits, hdu, tkeys, sizeof(tkeys)/sizeof(Fitkey), EOH_SEEK,
		   &key) == KEY_FOUND) {
/*
 * Install the keyword value in the appropriate slot.
 */
      if(key.extn > 0 && key.extn <= ahdu->tfields) {
	field = &ahdu->fields[key.extn-1];
	switch(key.keyid) {
	case AHDU_TBCOL:
	  field->tbcol = KEYINT(key);
	  break;
	case AHDU_TSCAL:
	  field->tscal = KEYDBL(key);
	  break;
	case AHDU_TZERO:
	  field->tzero = KEYDBL(key);
	  break;
	case AHDU_TFORM:
	  field->tform = fitsstr(KEYSTR(key));
	  get_format(KEYSTR(key), field);
	  break;
	case AHDU_TNULL:
	  field->tnull = fitsstr(KEYSTR(key));
	  break;
	case AHDU_TTYPE:
	  field->ttype = fitsstr(KEYSTR(key));
	  break;
	case AHDU_TUNIT:
	  field->tunit = fitsstr(KEYSTR(key));
	  break;
	};
      };
    };
  };
/*
 * Check that the mandatory field-descriptor keywords.
 */
  for(i=0; i<ahdu->tfields; i++) {
    field = &ahdu->fields[i];
    if(field->tbcol==0) {
      fprintf(stderr, "Missing TBCOL%d keyword\n", i+1);
      return 1;
    };
    if(field->type==DAT_NON || field->width==0 || field->tform==NULL) {
      fprintf(stderr, "Missing TFORM%d keyword\n", i+1);
      return 1;
    };
  };
  return 0;
}

/*.......................................................................
 * Return the size of an ASCII-TABLE HDU descriptor.
 *
 * Output:
 *  return  size_t  sizeof(Ahdu).
 */
static SIZEFN(siz_ahdu)
{
  return sizeof(Ahdu);
}

/*.......................................................................
 * Return the number of the ASCII-table-column that has a specified name.
 *
 * Input:
 *  thdu    thdu *  The generic table HDU descriptor.
 *  ttype   char *  The name of the column to be sought - trailing
 *                  spaces are ignored.
 *  fixlen  int     If > 0, then this defines the maximum number of characters
 *                  to be compared. This makes it possible to search using
 *                  prefixes.
 * Output:
 *  return   int    1-relative column number, or 0 if not found.
 */
static COL_FINDFN(acol_find)
{
  Ahdu *ahdu;  /* Ahdu version of thdu */
  Afield *af;  /* Array of ASCII field descriptors. */
  int field;   /* Field number being checked */
/*
 * Sanity check the HDU type.
 */
  if(thdu==NULL || thdu->type != F_TABLE) {
    fprintf(stderr, "acol_find: Bad HDU type received\n");
    return 0;
  };
/*
 * Get the descriptor and the array of field descriptors.
 */
  ahdu = (Ahdu *) thdu;
  af = ahdu->fields;
/*
 * Search out the column by the field names recorded in the field descriptors.
 */
  for(field=0; field<ahdu->tfields; field++) {
    if(matchstr(af[field].ttype, ttype, fixlen))
      return field+1;
  };
  return 0;
}

/*.......................................................................
 * Return the data-type of a given column.
 * NB. this routine should not be called directly - use the generic
 * col_type() function instead.
 *
 * Input:
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 * Output:
 *  return  Fittype    The data-type of the column. On error DAT_NON==0 is
 *                     returned.
 */
static COL_TYPEFN(acol_type)
{
  Afield *field = &((Ahdu *) thdu)->fields[icol-1];
/*
 * Get the data-type from the field descriptor.
 */
  return field->type;
}

/*.......................................................................
 * Return the dimension of a given table field.
 * NB. this routine should not be called directly - use the generic
 * col_dim() function instead.
 *
 * Input:
 *  fits       Fits *  The FITS file descriptor.
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 *  irow        int    The 1-relative row number. This is ignored.
 * Output:
 *  return      int    The number of elements across the specified column.
 *                     For ASCII tables this is always 1, except for
 *                     string fields where it is the number of characters
 *                     in the string.
 */
static COL_DIMFN(acol_dim)
{
  Afield *field = &((Ahdu *) thdu)->fields[icol-1];
/*
 * ASCII tables only support scalar columns.
 */
  return (field->type==DAT_CHR) ? field->width : 1;
}

/*.......................................................................
 * Decompose a format string to ascertain the field width and data-type.
 *
 * Input:
 *  tform     char *  The tform keyword - ignored if NULL.
 * Input/Output:
 *  field   Afield *  The column descriptor. 'tform' will be parsed
 *                    and the results placed in the 'type','width' and
 *                    'ndec' members of 'field'. On error these are
 *                    left unchanged (at their initialized values).
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int get_format(char *tform, Afield *field)
{
  Fittype type;   /* Temporary storage of the determined data-type */
  int width=0;    /* Temporary storage of the determined field width */
  int ndec=0;     /* Temporary storage of the number of decimal places */
  char *cptr;     /* Pointer into tformat */
  char *endp;     /* Pointer into tformat */
  int ierr=0;     /* Error status */
  if(field) {
/*
 * Determine the data-type from the first character of the tform keyword.
 */
    switch(*tform) {
    case 'A':
      type = DAT_CHR;
      break;
    case 'I':
      type = DAT_LNG;
      break;
    case 'F':
      type = DAT_FLT;
      break;
    case 'E':
      type = DAT_FLT;
      break;
    case 'D':
      type = DAT_DBL;
      break;
    default:
      type = DAT_NON;
      ierr = 1;
    };
/*
 * Determine the field width.
 */
    if(!ierr) {
      width = strtol(&tform[1], &cptr, 10);
      ierr = width==0;
    };
/*
 * Determine the number of decimal places.
 */
    if(!ierr && *cptr=='.') {
      ndec = strtol(++cptr, &endp, 10);
      ierr = cptr==endp;
    };
/*
 * If no error occurred record the details in the field descriptor.
 */
    if(ierr)
      fprintf(stderr, "get_format: Bad format: TFORM=%s\n", tform);
    else {
      field->type = type;
      field->width = width;
      field->ndec = ndec;
      field->form = *tform;
    };
  };
  return ierr;
}

/*.......................................................................
 * Return the name of a given column.
 * This function should be called via col_name(). Do not call directly
 * since no sanity checking is performed here in.
 *
 * Input:
 *  thdu     Thdu *  The generic table descriptor.
 *  icol      int    The number of the column to be looked up.
 * Output:
 *  return   char *  The column name or NULL if not known or on error.
 */
static COL_NAMEFN(acol_name)
{
  return ((Ahdu *) thdu)->fields[icol-1].ttype;
}

/*.......................................................................
 * Return values from a single table entry in a given column.
 *
 * Input:
 *  fits       Fits    The FITS file descriptor.
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 *  irow        int    The 1-relative row number.
 *  type    Fittype    The data-type of the 'data[]' array.  Implicit
 *                     conversion to this type will be performed if
 *                     meaningful.
 *  doscale     int    If true apply offset and scale factors.
 *  flags  Fitsflag *  If not NULL, then this should point to an array
 *                     of nobj elements to be used to record which elements
 *                     of data[] are blanked 1=blanked, 0=OK.
 *  first       int    This item is ignored since ASCII table column entries
 *                     are scalar.
 *  ndata       int    The max number of elements to return from the field.
 *                     This is ignored for all types except DAT_CHR where
 *                     it should be the number of char elements available
 *                     in the data[] array.
 * Output:
 *  data       void *  Supply an array of the type encoded in 'type' and
 *                     having at least 'ndata' elements.
 *  return     long    Number of objects of type 'type' read.
 *                     This is 0 on error and 1 for all types but DAT_CHR.
 */
static COL_VALFN(acol_value)
{
  Ahdu *ahdu = (Ahdu *) thdu;
  Afield *field;    /* The field descriptor */
  long offset;      /* 8-bit byte offset into data segment of HDU */
  char *tmpbuf=NULL;/* Pointer to malloc'd buffer to be freed */
  char *cvalue;     /* String holding the field value */
  char *cptr;       /* Pointer into cvalue string */
  int isnull=0;     /* True if the field is null */
  int ierr=0;       /* Error status */
/*
 * Create a union of all fits types.
 */
  static union {
    char c;
    int i;
    long l;
    float f;
    double d;
  } u;
/*
 * Get the field descriptor for column 'icol'.
 */
  field = &ahdu->fields[icol-1];
/*
 * Sanity check 'ndata'. For char-type fields, it must be constrained
 * to the field width. For other types ndata=1.
 */
  if(type==DAT_CHR) {
    if(ndata > field->width)
      ndata = field->width;
  } else {
    ndata = 1;
  };
/*
 * Get a buffer of at least the field width plus room for a
 * '\0' terminator. (We have to be able to read the whole field
 * even for DAT_CHR, so that TNULL can be checked against the field.
 */
  if(field->width+1 < MAX_ITEM_WIDTH)
    cvalue = item_buff;
  else
    cvalue = tmpbuf = calloc(field->width+1, sizeof(char));
/*
 * Buffer space unavailable?
 */
  if(cvalue==NULL) {
    fprintf(stderr,
        "acol_value: Insufficient memory to read a table field of %d bytes.\n",
	 field->width);
    return 0L;
  };
/*
 * Read the requested table field.
 */
  offset = (long) ahdu->dims[0] * (irow-1) + (field->tbcol-1);
  if(get_data(fits, (Hdu *)thdu, offset, DAT_CHR, 0L, field->width, DAT_CHR,
	      0.0, 1.0, NULL, NULL, NONULL, cvalue))
    return acol_return(tmpbuf, 0L);
/*
 * Terminate the string that was read.
 */
  cvalue[field->width] = '\0';
/*
 * Check for TNULL value.
 */
  if(field->tnull && matchstr(cvalue, field->tnull, 0)) {
    isnull = 1;
    if(flags) flags[0] = 1;
  };
/*
 * Read the value from the field.
 */
  switch(field->type) {
  case DAT_INT:
    u.i = isnull ? 0:strtol(cvalue, &cptr, 10);
    ierr = (cptr==cvalue && !isnull);
    break;
  case DAT_LNG:
    u.l = isnull ? 0L:strtol(cvalue, &cptr, 10);
    ierr = (cptr==cvalue && !isnull);
    break;
  case DAT_FLT:
    u.f = isnull ? 0.0f : (float) strtod(cvalue, &cptr);
    ierr = (cptr==cvalue && !isnull);
    break;
  case DAT_DBL:  /* FORTRAN 'D' format */
/*
 * Convert the 'D' exponent to an 'E' format so that C can read it.
 */
    cptr = strchr(cvalue, 'D');
    if(cptr) *cptr = 'E';
    u.d = isnull ? 0.0 : strtod(cvalue, &cptr);
    ierr = (cptr==cvalue && !isnull);
    break;
  case DAT_CHR:
    break;
  default:
    fprintf(stderr, "acol_value: Unhandled field data-type\n");
    ierr = 1;
  };
/*
 * Apply field zero offset and scale factors and convert to the
 * type of 'data'.
 */
  if(type==DAT_CHR)
    memcpy(data, cvalue, ndata);
  else
    ierr = ierr || typeconv(ndata, field->type, &u, doscale ? field->tzero:0.0,
			    doscale ? field->tscal:1.0, type, data);
/*
 * Error reading field?
 */
  if(ierr) {
    fprintf(stderr,
	    "Error reading field=\'%s\' at row %d of column %d of table %s\n",
	    cvalue, irow, icol, thdu->extname ? thdu->extname : "(no name)");
    return acol_return(tmpbuf, 0L);
  };
  return acol_return(tmpbuf, ndata);
}

/*.......................................................................
 * Private cleanup-on-return function of acol_set() and acol_value(),
 * used to de-allocate the optional extended I/O buffer and return a
 * given return code.
 *
 * Input:
 *  tmpbuf   char *   The buffer to be free'd (if not NULL).
 *  icode    long     The required return code.
 * Output:
 *  return   long     icode.
 */
static long acol_return(char *tmpbuf, long iret)
{
  if(tmpbuf) free(tmpbuf);
  return iret;
}

/*.......................................................................
 * Set values in a single ASCII-table column entry.
 *
 * Input:
 *  fits       Fits    The FITS file descriptor.
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 *  irow        int    The 1-relative row number.
 *  type    Fittype    The data-type of the 'data[]' array. Implicit
 *                     conversion from this type will be performed if
 *                     meaningful.
 *  doscale     int    If true remove offset and scale factors.
 *  flags Fitsflag *   If not NULL, then this should point to an array
 *                     of nobj elements specifying which elements
 *                     of data[] are to be blanked 1=blanked, 0=OK.
 *  first       int    This item is ignored since ASCII table column entries
 *                     are scalar.
 *  ndata       int    The max number of elements to write to the field.
 *                     This is ignored for all types except DAT_CHR where
 *                     it should be the number of char elements available
 *                     in the data[] array.
 *  data       void *  The array to be written. This must have type 'type'
 *                     and have at least 'ndata' elements.
 * Output:
 *  return     long    Number of objects of type 'type' written.
 *                     This is 0 on error and 1 for all types but DAT_CHR.
 */
static COL_SETFN(acol_set)
{
  const int safety=10;  /* Extra space in work-buffer for sprintf() overflows */
  Ahdu *ahdu = (Ahdu *) thdu;
  Afield *field;    /* The field descriptor */
  long offset;      /* 8-bit byte offset into data segment of HDU */
  long nobj;        /* The actual number of objects to write */
  char *tmpbuf=NULL;/* Pointer to malloc'd buffer to be freed */
  char *cvalue;     /* String holding the field value */
  char *cptr;       /* Pointer into output buffer */
  int i;
/*
 * Create a union of all fits types.
 */
  static union {
    char c;
    int i;
    long l;
    float f;
    double d;
  } u;
/*
 * Get the field descriptor for column 'icol'.
 */
  field = &ahdu->fields[icol-1];
/*
 * Sanity check 'ndata'. For char-type fields, it must be constrained
 * to the field width. For other types ndata=1.
 */
  if(type==DAT_CHR) {
    if(ndata > field->width)
      ndata = field->width;
  } else {
    ndata = 1;
  };
/*
 * Except when TNULL must be substituted, nobj equals ndata.
 */
  nobj = ndata;
/*
 * Get a buffer of at least the field width plus room for a
 * '\0' terminator plus room for sprintf() overflows.
 */
  if(field->width+safety < MAX_ITEM_WIDTH)
    cvalue = item_buff;
  else
    cvalue = tmpbuf = calloc(field->width+safety, sizeof(char));
/*
 * Buffer space unavailable?
 */
  if(cvalue==NULL) {
    fprintf(stderr,
        "acol_set: Insufficient memory to write a table field of %d bytes.\n",
	 field->width);
    return 0L;
  };
/*
 * Substitute field->tnull for data[] if the field is blanked.
 */
  if(flags && flags[0] && field->tnull) {
    nobj = field->width;
    data = field->tnull;
    type = DAT_CHR;
  };
/*
 * Convert from the type given to the declared type of the field.
 */
  if(type != DAT_CHR) {
    if(typeconv(nobj, type, data, doscale ? -field->tzero:0.0,
		doscale ? 1.0/field->tscal:1.0, field->type, &u))
      return acol_return(tmpbuf, 0L);
  };
/*
 * Write the field value into cvalue[].
 */
  switch(field->type) {
  case DAT_INT:
    sprintf(cvalue, "%*d", (int)field->width, u.i);
    break;
  case DAT_LNG:
    sprintf(cvalue, "%*ld", (int)field->width, u.l);
    break;
  case DAT_FLT:
    if(field->form=='E')
      sprintf(cvalue, "%*.*E", (int)field->width, field->ndec, u.f);
    else
      sprintf(cvalue, "%*.*f", (int)field->width, field->ndec, u.f);
    break;
  case DAT_DBL:  /* FORTRAN 'D' format */
    sprintf(cvalue, "%*.*E", (int)field->width, field->ndec, u.d);
/*
 * Convert the 'E' exponent to a 'D'.
 */
    cptr = strchr(cvalue, 'E');
    if(cptr) *cptr = 'D';
    break;
  case DAT_CHR:
/*
 * Copy and terminate the string.
 */
    memcpy(cvalue, data, nobj);
    cvalue[nobj] = '\0';
    break;
  default:
    fprintf(stderr, "acol_set: Unhandled ASCII table field-data-type.\n");
    return 0L;
  };
/*
 * Space pad the output field.
 */
  cvalue[field->width]='\0';
  for(i=strlen(cvalue); i<field->width; i++)
    cvalue[i]=' ';
  cvalue[i] = '\0';
/*
 * Write the requested table field.
 */
  offset = (long) ahdu->dims[0] * (irow-1) + (field->tbcol-1);
  if(put_data(fits, (Hdu *)thdu, offset, DAT_CHR, 0L, field->width, DAT_CHR,
	      0.0,1.0,NULL, NULL, NONULL, cvalue))
    return acol_return(tmpbuf, 0L);
  return acol_return(tmpbuf, ndata);
}

/*.......................................................................
 * Create a new ASCII-table extension HDU descriptor.
 *
 * Input:
 *  width    int   The number of characters across one row of the table.
 *  nrow     int   The number of rows in the table.
 *  extname char * A name to refer to this extension by.
 *                 A copy of the string will be made.
 *  extver   int   The version number to assign to the extension. If
 *                 0 - then when the HDU is added to a FITS file, it
 *                 will be assigned the next highest version number for
 *                 name 'extname' and level 'extlevel'.
 *  extlevel int   The hierachical level of this extension. Send 1 or
 *                 0 if not relevant.
 *  tfields  int   The number of fields per row.
 * Output:
 *  return   int   0 - OK.
 *                 1 - Error.
 */ 
Hdu *new_asctab(int width, int nrow, char *extname, int extver,	int extlevel,
		int tfields)
{
  Hdu *hdu;
  Ahdu *ahdu;
  int dims[2];
/*
 * Allocate and default initialize a table HDU.
 */
  hdu = new_Hdu(F_TABLE);
  if(hdu==NULL)
    return hdu;
/*
 * Get the derived version of hdu.
 */
  ahdu = (Ahdu *) hdu;
/*
 * Initialize the HDU base-class fields.
 */
  dims[0] = width;
  dims[1] = nrow;
  hdu = ini_Hdu(hdu, B_CHAR, dims, 2, 0, 0, 1, extname, extver, extlevel,
		0, 0);
/*
 * Assign the table base-class fields.
 */
  if(tfields < 1) {
    fprintf(stderr, "new_asctab: Illegal tfields=%d\n", tfields);
    return del_Hdu(hdu);
  };
  ahdu->tfields = tfields;
/*
 * Allocate the array of field descriptors.
 */
  ahdu->fields = new_Afields(ahdu);
  if(ahdu->fields==NULL && ahdu->tfields>0)
    return del_Hdu(hdu);
  return hdu;
}

/*.......................................................................
 * Define the structure of a new table column.
 *
 * Input:
 *  hdu      Hdu *  The descriptor of the new table HDU.
 *  icol     int    The 1-relative field number in the table.
 *  tbcol    int    The start character of the new column in the table.
 *  tscal double    The scale factor for the new field.
 *  tzero double    The zero offset for the new field.
 *  tform   char *  The mandatory FORTRAN-77 format for the field.
 *                  The only valid formats are:
 *                     Aw  Character array of size w.
 *                     Iw  Integer in a field width of w characters.
 *                   Fw.d  Float in a field width of w chars with
 *                         d decimal places displayed.
 *                   Ew.d  As Fw.d but using E exponential notation.
 *                   Dw.d  As Fw.d but for double-precision and displayed
 *                         using FORTRAN D exponential notation.
 *  tnull  char *   An optional string used to denote NULL values in this
 *                  field. A string copy - implicitly filled to the field
 *                  width with blanks will be made.
 *  ttype  char *   An optional name for this field.
 *  tunit  char *   An optional name for the units of measurement of this
 *                  field.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Error.
 */
int setafield(Hdu *hdu, int icol, int tbcol, double tscal, double tzero,
	      char *tform, char *tnull, char *ttype, char *tunit)
{
  Ahdu *ahdu = (Ahdu *) hdu;
  Afield *field;
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL || hdu->type != F_TABLE) {
    fprintf(stderr, "setafield: Inappropriate HDU descriptor received\n");
    return 1;
  };
  if(hdu->complete) {
    fprintf(stderr,
	  "setafield: Illegal attempt to change established HDU structure\n");
    return 1;
  };
  if(tform==NULL) {
    fprintf(stderr, "setafield: Missing TFORM argument\n");
    return 1;
  };
  if(icol<1 || icol>ahdu->tfields) {
    fprintf(stderr, "setafield: Column index out of range\n");
    return 1;
  };
  if(tbcol < 1 || tbcol > ahdu->dims[0]) {
    fprintf(stderr, "setafield: tbcol argument out of range\n");
    return 1;
  };
  if(tscal==0.0) {
    fprintf(stderr, "setafield: Error: tscal=0\n");
    return 1;
  };
/*
 * Get the field descriptor.
 */
  field = &ahdu->fields[icol-1];
/*
 * Decompose the tform keyword.
 */
  if(get_format(tform, field))
    return 1;
/*
 * Assign the rest of the members.
 */
  field->tbcol = tbcol;
  field->tscal = tscal;
  field->tzero = tzero;
  if(!field->tform && tform)
    field->tform = fitsstr(tform);
  if(!field->tnull && tnull)
    field->tnull = fitsstr(tnull);
  if(!field->ttype && ttype)
    field->ttype = fitsstr(ttype);
  if(!field->tunit && tunit)
    field->tunit = fitsstr(tunit);
  return 0;
}

/*.......................................................................
 * Write the header lines for the derived parts of an ASCII table
 * extension HDU to a FITS file.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file to which the HDU is
 *                  being added.
 *  hdu      Hdu *  An initialized HDU descriptor for the new HDU.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static ADDFN(add_ahdu)
{
  Ahdu *ahdu = (Ahdu *) hdu;
  Afield *field;
  int i;
/*
 * Check that all table fields have been initialized.
 */
  field = ahdu->fields;
  for(i=0; i<ahdu->tfields; i++, field++) {
    if(field->type == DAT_NON) {
      fprintf(stderr, "add_ahdu: HDU incomplete: Please use setafield()"
	      " to describe table field %d\n", i+1);
      return 1;
    };
  };
/*
 * Write the Mandatory PCOUNT, GCOUNT and TFIELDS keywords.
 */
  if(wintkey(fits, hdu, "PCOUNT", 0, hdu->pcount, "Parameter count") ||
     wintkey(fits, hdu, "GCOUNT", 0, hdu->gcount, "Group count") ||
     wintkey(fits, hdu, "TFIELDS",0, ahdu->tfields, "Number of table fields"))
    return 1;
/*
 * Write extension name and version keywords where given.
 */
  if(w_extkeys(fits, hdu))
    return 1;
/*
 * Write the field description keywords where given.
 */
  if(ahdu->fields) {
    for(i=0; i<ahdu->tfields; i++) {
      field = &ahdu->fields[i];
      if(field->ttype && wstrkey(fits, hdu, "TTYPE", i+1, field->ttype,
				 "Name of this table field"))
	return 1;
      if(field->tunit && wstrkey(fits, hdu, "TUNIT", i+1, field->tunit,
				 "Unit of measurement of this table field"))
	return 1;
      if(field->tnull && wstrkey(fits, hdu, "TNULL", i+1, field->tnull,
				 "Value used to indicate a NULL item"))
	return 1;
      if(field->tform && wstrkey(fits, hdu, "TFORM", i+1, field->tform,
				 "Format of table field"))
	return 1;
      if(wintkey(fits, hdu, "TBCOL", i+1, field->tbcol,
		 "Start character in row"))
	return 1;
      if(field->tscal!=1.0 && wfltkey(fits, hdu, "TSCAL", i+1, field->tscal,
		 "Scale factor applied to items in this field"))
	return 1;
      if(field->tzero!=0.0 && wfltkey(fits, hdu, "TZERO", i+1, field->tzero,
		 "Zero offset applied to items in this field"))
	return 1;
    };
  };
  return 0;
}

/*.......................................................................
 * Allocate ahdu->tfields ASCII-table descriptors.
 *
 * Input:
 *  ahdu     Ahdu *  The descriptor of the ASCII table HDU in which to
 *                   allocate the array.
 * Output:
 *  return Afield *  The allocated array of descriptors.
 */
static Afield *new_Afields(Ahdu *ahdu)
{
  Afield *field;
  int i;
/*
 * Sanity checks.
 */
  if(ahdu==NULL) {
    fprintf(stderr, "new_Afields: Intercepted NULL Ahdu descriptor\n");
    return NULL;
  };
  if(ahdu->fields) {
    fprintf(stderr, "new_Afields: ahdu->fields allready allocated\n");
    return NULL;
  };
/*
 * Allocate an array of table field descriptors.
 */
  if(ahdu->tfields>0) {
    ahdu->fields = (Afield *) malloc(sizeof(Afield) * ahdu->tfields);
    if(ahdu->fields == NULL) {
      fprintf(stderr,
       "new_Afields: Insufficient memory for ASCII-table field descriptors\n");
      return NULL;
    };
    for(i=0; i<ahdu->tfields; i++) { /* Initialize members */
      field = &ahdu->fields[i];
      field->type  = DAT_NON;
      field->tbcol = 0;
      field->tscal = 1.0;
      field->tzero = 0.0;
      field->width = 0;
      field->ndec  = 0;
      field->form  = ' ';
      field->tform = NULL;
      field->tnull = NULL;
      field->ttype = NULL;
      field->tunit = NULL;
    };
  };
  return ahdu->fields;
}

/*.......................................................................
 * Copy the derived parts of a Ahdu descriptor into another Ahdu
 * descriptor. This function should be called only via copy_Hdu().
 *
 * Input:
 *  hdu    Hdu * The original hdu to be copied.
 * Output:
 *  return Hdu * The new copied version of 'hdu', or NULL on error.
 */
static COPYFN(cop_ahdu)
{
  Ahdu *old = (Ahdu *) hdu;
  Hdu *new;   /* The new HDU descriptor */
  Afield *f;  /* A field descriptor from the original HDU descriptor */
  int icol;   /* Column number in table */
/*
 * Create the new ASCII-table descriptor.
 */
  new = new_asctab(old->dims[0], old->dims[1], old->extname, 0, old->extlevel,
		   old->tfields);
  if(new==NULL)
    return new;
/*
 * Copy the field descriptors.
 */
  f = old->fields;
  for(icol=1; icol<=old->tfields; icol++, f++) {
    if(setafield(new, icol, f->tbcol, f->tscal, f->tzero, f->tform, f->tnull,
		 f->ttype, f->tunit))
      return del_Hdu(new);
  };
  return 0;
}

/*.......................................................................
 * Complete the data section of an ASCII table HDU.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file in which the HDU resides.
 *  hdu      Hdu *  The descriptor of the HDU to be completed.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static ENDFN(end_ahdu)
{
  int saveline; /* The current header line number */
  int waserr=0; /* Error status */
/*
 * Now that the number of rows is known, re-write the NAXIS2
 * header keyword.
 */
  saveline = new_hline(hdu, 4);
  waserr=wintkey(fits, hdu, "NAXIS", 2, hdu->dims[1], "Number of table rows.");
  new_hline(hdu, saveline);
  return waserr!=0;
}

