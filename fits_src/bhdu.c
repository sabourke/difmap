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
#include "bhdu.h"
#include "fitkey.h"

static int get_format(char *tform, Bfield *field);
static Bfield *new_Bfields(Bhdu *bhdu);

static COL_FINDFN(bcol_find);
static COL_TYPEFN(bcol_type);
static COL_DIMFN(bcol_dim);
static COL_VALFN(bcol_value);
static COL_NAMEFN(bcol_name);
static COL_SETFN(bcol_set);

Tabfn btabfn={bcol_value, bcol_find, bcol_type, bcol_dim, bcol_name, bcol_set};

static GETFN(get_bhdu);
static NEWFN(new_bhdu);
static DELFN(del_bhdu);
static SIZEFN(siz_bhdu);
static ADDFN(add_bhdu);
static COPYFN(cop_bhdu);
static ENDFN(end_bhdu);

Hdutab bhdufns={"BINTABLE", "A3DTABLE",	get_bhdu, new_bhdu, del_bhdu,
		siz_bhdu, add_bhdu, cop_bhdu, end_bhdu};

/*.......................................................................
 * Clear the derived part of an BINARY-table header descriptor to the
 * point at which it can be safely sent to del_Hdu().
 *
 * Input:
 *  hdu     Hdu *  The base-class pointer to a Bhdu descriptor.
 * Output:
 *  return  Hdu *  'hdu', or NULL if an error occurs and hdu is deleted.
 */
static NEWFN(new_bhdu)
{
  Bhdu *bhdu = (Bhdu *) hdu;
  (void) new_table(hdu);
  bhdu->fields = NULL;
  bhdu->theap = 0;
  bhdu->heap_nxt = 0;
  return hdu;
}

/*.......................................................................
 * Delete the derived parts of an BINARY-table HDU descriptor.
 *
 * Input:
 *  hdu      Hdu *  The base-class pointer to the Bhdu descriptor.
 */
static DELFN(del_bhdu)
{
  Bhdu *bhdu = (Bhdu *) hdu;
  Bfield *field;
  int i;
/*
 * Delete the generic table members.
 */
  del_table(hdu);
/*
 * Clean-up the array of table field descriptors.
 */
  if(bhdu->fields) {
    for(i=0; i<bhdu->tfields; i++) {
      field = &bhdu->fields[i];
      if(field->tform)
	free(field->tform);
      if(field->ttype)
	free(field->ttype);
      if(field->tunit)
	free(field->tunit);
      if(field->tdisp)
	free(field->tdisp);
      if(field->tdim)
	free(field->tdim);
    };
    free(bhdu->fields);
  };
  return;
}

/*.......................................................................
 * Read the BINARY-table parts of a FITS header and record
 * the table description in an BINARY-table-HDU descriptor.
 *
 * Input:
 *  fits  Fits *  The descriptor of the FITS file.
 * Input/Output:
 *  hdu    Hdu *  The base-class descriptor to a Bhdu descriptor.
 * Output:
 *  return int    0 - OK.
 */
static GETFN(get_bhdu)
{
  Bhdu *bhdu = (Bhdu *) hdu;
  Bfield *field;     /* A field descriptor */
  Fitkey key;        /* Keyword-value pair temporary storage */
  size_t tbcol;      /* Bytewize column in table */
  int i;
  /* Table field descriptor keywords */
  enum {BHDU_TFORM, BHDU_TSCAL, BHDU_TZERO, BHDU_TNULL, BHDU_TTYPE,
        BHDU_TUNIT, BHDU_TDISP, BHDU_TDIM};
  Fitkey tkeys[]={
    {"TFORM", 0, BHDU_TFORM, DAT_STR, NULL, NULL},
    {"TSCAL", 0, BHDU_TSCAL, DAT_DBL, NULL, NULL},
    {"TZERO", 0, BHDU_TZERO, DAT_DBL, NULL, NULL},
    {"TNULL", 0, BHDU_TNULL, DAT_INT, NULL, NULL},
    {"TTYPE", 0, BHDU_TTYPE, DAT_STR, NULL, NULL},
    {"TUNIT", 0, BHDU_TUNIT, DAT_STR, NULL, NULL},
    {"TDISP", 0, BHDU_TDISP, DAT_STR, NULL, NULL},
    {"TDIM",  0, BHDU_TDIM,  DAT_STR, NULL, NULL},
  };
/*
 * Check the input Hdu type.
 */
  if(hdu->type != F_BINTAB) {
    fprintf(stderr, "get_bhdu: Incompatible HDU descriptor received\n");
    return 1;
  };
/*
 * Check that NAXIS has the expected value.
 */
  if(hdu->naxis != 2) {
    fprintf(stderr, "Invalid NAXIS value (%d) != 2 in a binary table\n",
	    hdu->naxis);
    return 1;
  };
/*
 * The next keywords should be PCOUNT and GCOUNT - get_Hdu will have
 * already read there values so there is no need to assign them again here.
 */
  if(get_key(fits, hdu, "PCOUNT", DAT_INT, NO_SEEK, &key))
    fprintf(stderr, "Missing PCOUNT keyword in FITS table header\n");
  if(get_key(fits, hdu, "GCOUNT", DAT_INT, NO_SEEK, &key))
    fprintf(stderr, "Missing GCOUNT keyword in FITS table header\n");
/*
 * Check the sanity of the PCOUNT and GCOUNT values.
 */
  if(bhdu->gcount != 1 || bhdu->pcount < 0) {
    fprintf(stderr, "Illegal values of PCOUNT=%d GCOUNT=%d in binary table\n",
	    bhdu->pcount, bhdu->gcount);
    return 1;
  };
/*
 * Get the TFIELDS keyword (number of table fields per row).
 */
  if(get_key(fits, hdu, "TFIELDS", DAT_INT, NO_SEEK, &key)) {
    fprintf(stderr, "Missing TFIELDS keyword in FITS table header\n");
  } else {
    bhdu->tfields = KEYINT(key);
  };
/*
 * Allocate an array of table field descriptors.
 */
  if(new_Bfields(bhdu)==NULL)
    return 1;
/*
 * Get the table-field descriptor keywords. The TFORM keyword
 * is the only mandatory keyword.
 */
  if(bhdu->tfields > 0) {
    while(next_key(fits, hdu, tkeys, sizeof(tkeys)/sizeof(Fitkey), EOH_SEEK,
		   &key) == KEY_FOUND) {
/*
 * Install the keyword value in the appropriate slot.
 */
      if(key.extn > 0 && key.extn <= bhdu->tfields) {
	field = &bhdu->fields[key.extn-1];
	switch(key.keyid) {
	case BHDU_TSCAL:
	  field->tscal = KEYDBL(key);
	  break;
	case BHDU_TZERO:
	  field->tzero = KEYDBL(key);
	  break;
	case BHDU_TFORM:
	  field->tform = fitsstr(KEYSTR(key));
	  get_format(KEYSTR(key), field);
	  break;
	case BHDU_TNULL:
	  field->tnull = KEYINT(key);
	  break;
	case BHDU_TTYPE:
	  field->ttype = fitsstr(KEYSTR(key));
	  break;
	case BHDU_TUNIT:
	  field->tunit = fitsstr(KEYSTR(key));
	  break;
	case BHDU_TDISP:
	  field->tdisp = fitsstr(KEYSTR(key));
	  break;
	case BHDU_TDIM:
	  field->tdim = fitsstr(KEYSTR(key));
	  break;
	};
      };
    };
  };
/*
 * Check that the mandatory field-descriptor keywords.
 */
  for(i=0; i<bhdu->tfields; i++) {
    field = &bhdu->fields[i];
    if(field->type==DAT_NON || field->tform==NULL) {
      fprintf(stderr, "Missing TFORM%d keyword\n", i+1);
      return 1;
    };
  };
/*
 * Determine field start columns from the widths of each field.
 */
  tbcol = 1;
  field = bhdu->fields;
  for(i=0; i<bhdu->tfields; i++,field++) {
    field->tbcol = tbcol;
    tbcol += field->width * (field->isvar ? 8 : typesize(field->type));
  };
/*
 * Search for the optional THEAP keyword.
 */
  if(get_key(fits, hdu, "THEAP", DAT_INT, LOOP_SEEK, &key))
    bhdu->theap = bhdu->dims[0] * bhdu->dims[1];
  else
    bhdu->theap = KEYINT(key);
  return 0;
}

/*.......................................................................
 * Return the size of a BINARY-TABLE HDU descriptor.
 *
 * Output:
 *  return  size_t  sizeof(Bhdu).
 */
static SIZEFN(siz_bhdu)
{
  return sizeof(Bhdu);
}

/*.......................................................................
 * Return the number of the BINARY-table-column that has a specified name.
 *
 * Input:
 *  thdu    Thdu *  The generic table HDU descriptor.
 *  ttype   char *  The name of the column to be sought - trailing spaces
 *                  are ignored in the comparisons.
 *  fixlen  int     If > 0, then this defines the maximum number of characters
 *                  to be compared. This makes it possible to search by
 *                  prefixes.
 * Output:
 *  return   int    1-relative column number, or 0 if not found.
 */
static COL_FINDFN(bcol_find)
{
  Bhdu *bhdu;  /* Bhdu version of thdu */
  Bfield *bf;  /* Array of BINARY field descriptors. */
  int field;   /* Field number being checked */
/*
 * Sanity check the HDU type.
 */
  if(thdu==NULL || thdu->type != F_BINTAB) {
    fprintf(stderr, "bcol_find: Bad HDU type received\n");
    return 0;
  };
/*
 * Get the descriptor and the array of field descriptors.
 */
  bhdu = (Bhdu *) thdu;
  bf = bhdu->fields;
/*
 * Compare each field name with the required string.
 */
  for(field=0; field<bhdu->tfields; field++) {
    if(matchstr(bf[field].ttype, ttype, fixlen))
      return field+1;
  };
  return 0;
}

/*.......................................................................
 * Decompose a format string to ascertain the field dimensions and
 * data-type.
 *
 * Input:
 *  tform     char *  The tform keyword - ignored if NULL.
 * Input/Output:
 *  field   Bfield *  The column descriptor. 'tform' will be parsed
 *                    and the results placed in the 'type' and
 *                    'width' members of 'field'. On error these are
 *                    left unchanged (at their initialized values).
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int get_format(char *tform, Bfield *field)
{
  Fittype type=DAT_NON; /* Temporary storage of the determined data-type */
  int width;      /* Temporary storage of the determined field width */
  int isvar=0;    /* True if variable length array */
  char *cptr;     /* Pointer into tformat */
  int ierr=0;     /* Error status */
  if(field) {
/*
 * Read the element count.
 */
    width = strtol(tform, &cptr, 10);
    ierr = tform==cptr;
/*
 * Is this a variable length array entry?
 */
    if(!ierr && *cptr=='P') {
      isvar = 1;
      cptr++;
    };
/*
 * Determine the data-type.
 */
    if(!ierr) {
      field->form = *cptr++;
      switch(field->form) {
      case 'L':
	type = DAT_LOG;
	break;
      case 'X':
	type = DAT_BIT;
	break;
      case 'B':
	type = DAT_BYT;
	break;
      case 'I':
	type = DAT_SHT;
	break;
      case 'J':
	type = (CHAR_BIT*sizeof(int)>=32) ? DAT_INT : DAT_LNG;
	break;
      case 'A':
	type = DAT_CHR;
	break;
      case 'E':
	type = DAT_FLT;
	break;
      case 'D':
	type = DAT_DBL;
	break;
      case 'C':
	type = DAT_SCMP;
	break;
      case 'M':
	type = DAT_DCMP;
	break;
      default:
	type = DAT_NON;
	ierr = 1;
      };
    };
/*
 * If no error occurred record the details in the field descriptor.
 */
    if(ierr)
      fprintf(stderr, "get_format: Bad format: TFORM=%s\n", tform);
    else {
      field->type = type;
      field->width = width;
      field->isvar = isvar;
    };
  };
  return ierr;
}

/*.......................................................................
 * Return the data-type of a given column.
 * NB. this routine should not be called directly - use the generic
 * col_type() function instead.
 *
 * Input:
 *  fits       Fits *  The FITS file descriptor.
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 *  irow        int    The 1-relative row number. This is ignored except
 *                     when the column fields hold variable length arrays.
 * Output:
 *  return  Fittype    The data-type of the column. On error DAT_NON==0 is
 *                     returned.
 */
static COL_TYPEFN(bcol_type)
{
  Bhdu *bhdu = (Bhdu *) thdu;
  Bfield *field;
/*
 * Get the field descriptor.
 */
  field = &bhdu->fields[icol-1];
/*
 * Get the data-type from the field descriptor.
 */
  return field->type;
}

/*.......................................................................
 * Return the number of elements across a given column in a given row.
 * NB. this routine should not be called directly - use the generic
 * col_format() function instead.
 *
 * Input:
 *  fits       Fits *  The FITS file descriptor.
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 *  irow        int    The 1-relative row number. This is ignored except
 *                     when the column fields hold variable length arrays.
 * Output:
 *  return      int    The number of elements across the specified column
 *                     (and at the specified row for variable length arrays).
 *                     required. NB. This may be zero and is set to zero on
 *                     error or if irow or icol are out of bounds.
 */
static COL_DIMFN(bcol_dim)
{
  Bhdu *bhdu = (Bhdu *) thdu;
  Bfield *field;
  int ndata;
/*
 * Get the field descriptor.
 */
  field = &bhdu->fields[icol-1];
/*
 * Determine the number of elements in the field.
 */
  if(!field->isvar) {
    ndata = field->width;
  } else {
    long offset;   /* Offset of field wrt start of HDU data segment */
    long vdesc[2]; /* Variable length array descriptor */
/*
 * Check that the row index is within bounds.
 */
    if(irow<1 || irow>thdu->dims[1]) {
      ndata = 0;
    } else {
/*
 * Compute the offset to the field and read the variable length array
 * descriptor from that field.
 */
      offset = (long) bhdu->dims[0] * (irow-1) + (field->tbcol-1);
      if(get_data(fits, (Hdu *) thdu, offset, DAT_LNG, 0L, 2L,
		  DAT_LNG, 0.0, 1.0, NULL, NULL, NONULL, &vdesc))
	return 0;
      ndata = vdesc[0];                 /* Number of elements stored */
    };
  };
  return ndata;
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
static COL_NAMEFN(bcol_name)
{
  return ((Bhdu *) thdu)->fields[icol-1].ttype;
}

/*.......................................................................
 * Return values from a single table entry in a given column.
 *
 * Input:
 *  fits       Fits    The FITS file descriptor.
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 *  irow        int    The 1-relative row number.
 *  type    Fittype    The data-type of the 'data[]' array. Implicit
 *                     conversion to this type will be performed if
 *                     meaningful.
 *  doscale     int    If true apply offset and scale factors.
 *  flags  Fitsflag *  If not NULL, then this should point to an array
 *                     of nobj elements to be used to record which elements
 *                     of data[] are blanked 1=blanked, 0=OK.
 *  first       int    The 0-relative index of the first element to get.
 *  ndata       int    The max number of elements to return from the field.
 * Output:
 *  data       void *  Supply an array of the type encoded in 'type' and
 *                     having at least 'ndata' elements.
 *  return     long    Number of objects read. This will be less than
 *                     ndata if the field contains fewer than ndata objects.
 *                     On error 0 is returned.
 */
static COL_VALFN(bcol_value)
{
  Bhdu *bhdu = (Bhdu *) thdu;
  Bfield *field;    /* The field descriptor */
  long offset;      /* 8-bit byte offset into data segment of HDU */
  long nmax;        /* Number of elements available in the given field */
  long vdesc[2];    /* Variable array FITS field descriptor */
/*
 * Get the field descriptor for column 'icol'.
 */
  field = &bhdu->fields[icol-1];
/*
 * Compute the offset to the requested field.
 */
  offset = (long) bhdu->dims[0] * (irow-1) + (field->tbcol-1);
/*
 * If this is a variable length array - read the heap offset and
 * data count from the field.
 */
  if(field->isvar) {
    if(bhdu->state != HDU_INFILE) {
      fprintf(stderr,
     "wcolumn: You must call end_hdu() before using variable length arrays.\n");
      return 0;
    };
    if(get_data(fits, (Hdu *) thdu, offset, DAT_LNG, 0L, 2L,
		DAT_LNG, 0.0, 1.0, NULL, NULL, NONULL, &vdesc))
      return 0;
    nmax = vdesc[0];                 /* Number of elements stored */
    offset = bhdu->theap + vdesc[1]; /* Position in heap for read */
  } else {
    nmax = field->width;
  };
/*
 * Limit the number of objects to be read to the number of elements in
 * the field.
 */
  if(ndata+first > nmax)
    ndata = nmax-first;
/*
 * Read the requested table field.
 */
  if(get_data(fits, (Hdu *) thdu, offset, field->type, first, ndata,
	      type, doscale ? field->tzero:0.0, doscale ? field->tscal:1.0,
	      NULL, flags, field->tnull, data))
    return 0;
  else
    return ndata;
}

/*.......................................................................
 * Set values in a single table column entry. IMPORTANT: If you are
 * writing to a variable length array field then you must have pre-set
 * its dimension with a single previous call to
 * setdim(fits, thdu, icol, irow, dim). If you fail to do this then the
 * data will be written at a random position in the file and corrupt it.
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
 *  first       int    The 0-relative index of the first element to write.
 *  ndata       int    The max number of elements to write to the field.
 *  data       void *  The array to be written. This must have type 'type'
 *                     and have at least 'ndata' elements.
 * Output:
 *  return     long    Number of objects written. This will be less than
 *                     ndata if the field contains fewer than ndata objects,
 *                     or if type==DAT_CHR and data[] contains a string that
 *                     is terminated before the ndata'th character.
 *                     On error 0 is returned.
 */
static COL_SETFN(bcol_set)
{
  Bhdu *bhdu = (Bhdu *) thdu;
  Bfield *field;    /* The field descriptor */
  long offset;      /* 8-bit byte offset into data segment of HDU */
  long nmax;        /* Number of elements available in the given field */
  long vdesc[2];    /* Variable array FITS field descriptor */
  char *cptr;
  int i;
/*
 * Get the field descriptor for column 'icol'.
 */
  field = &bhdu->fields[icol-1];
/*
 * Compute the offset to the requested field.
 */
  offset = (long) bhdu->dims[0] * (irow-1) + (field->tbcol-1);
/*
 * If this is a variable length array - read the heap offset and
 * data count from the field.
 */
  if(field->isvar) {
    if(bhdu->state != HDU_INFILE) {
      fprintf(stderr,
     "rcolumn: You must call end_hdu() before using variable length arrays.\n");
      return 0;
    };
    if(get_data(fits, (Hdu *) thdu, offset, DAT_LNG, 0L, 2L,
		DAT_LNG, 0.0, 1.0, NULL, NULL, NONULL, &vdesc))
      return 0;
    nmax = vdesc[0];                 /* Number of elements stored */
    offset = bhdu->theap + vdesc[1]; /* Position in heap for read */
  } else {
    nmax = field->width;
  };
/*
 * Limit the number of objects to be read to the number of elements in
 * the field.
 */
  if(ndata+first > nmax)
    ndata = nmax-first;
/*
 * If a string is being written only write up to ndata or up to the
 * first '\0' terminator - whichever comes first.
 */
  if(type==DAT_CHR) {
    cptr = (char *) data;
    for(i=0; i<ndata && *cptr; i++, cptr++);
/*
 * If a '\0' character appears within the 'ndata' limit, include it in the
 * output to be written.
 */
    if(i<ndata)
      ndata = i+1;
  };
/*
 * Write the requested table field.
 */
  if(put_data(fits, (Hdu *) thdu, offset, field->type, first, ndata,
	      type, doscale ? -field->tzero:0.0, doscale ? 1.0/field->tscal:1.0,
	      NULL, flags, field->tnull, data))
    return 0L;
  else
    return ndata;
}

/*.......................................................................
 * Set the dimensions of a new variable length array table field. This
 * includes allocating space within the heap.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  thdu    Thdu *   The BINARY table descriptor.
 *  icol     int     The 1-relative number of the column.
 *  irow     int     The 1-relative row number.
 *  ndata    int     The number of elements to assign to the field.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int setdim(Fits *fits, Thdu *thdu, int icol, int irow, int ndata)
{
  Bhdu *bhdu = (Bhdu *) thdu;
  Bfield *field;    /* The field descriptor */
  long offset;      /* 8-bit byte offset into data segment of HDU */
  long vdesc[2];    /* Variable array FITS field descriptor */
  size_t size;      /* Size of an element in the variable length array */
  size_t nleft;     /* Number of remaining bytes of heap storage */
/*
 * Sanity check the arguments.
 */
  if(thdu==NULL) {
    fprintf(stderr, "setdim: NULL HDU descriptor received\n");
    return 1;
  };
  if(thdu->type != F_BINTAB) {
    fprintf(stderr, "setdim: Not a binary table!\n");
    return 1;
  };
  if(icol < 1 || icol > thdu->tfields) {
    fprintf(stderr, "setdim: Out-of-range column index received\n");
    return 1;
  };
  if(irow < 1 || irow > thdu->dims[1]) {
    fprintf(stderr, "setdim: Out-of-range row index received\n");
    return 1;
  };
/*
 * Get the field descriptor for column 'icol'.
 */
  field = &bhdu->fields[icol-1];
/*
 * Is it declared to hold variable length arrays?
 */
  if(!field->isvar) {
    fprintf(stderr, "setdim: Column %d does not hold variable length arrays\n",
	    icol);
    return 1;
  };
/*
 * Determine the size of type in FITS bytes.
 */
  size = typesize(field->type);
/*
 * How much room is left in the heap?
 */
  nleft = ((long) bhdu->dims[0]*bhdu->dims[1] + bhdu->pcount) -
          (bhdu->theap + bhdu->heap_nxt);
/*
 * Is there sufficient room in the heap?
 */
  if(nleft < size * ndata) {
    fprintf(stderr, "setdim: Insufficient room in heap for new array\n");
    return 1;
  };
/*
 * Prepare the variable array descriptor to be written and increment
 * the heap pointer.
 */
  vdesc[0] = ndata;
  vdesc[1] = bhdu->heap_nxt;
  bhdu->heap_nxt += ndata * size;
/*
 * Compute the offset to the requested field.
 */
  offset = (long) bhdu->dims[0] * (irow-1) + (field->tbcol-1);
/*
 * Write the descriptor field.
 */
  if(put_data(fits, (Hdu *) thdu, offset, DAT_LNG, 0L, 2L,
	      DAT_LNG, 0.0, 1.0, NULL, NULL, NONULL, &vdesc))
    return 1;
  return 0;
}

/*.......................................................................
 * Public function used to inquire whether a field contains a variable
 * length array.
 *
 * Input:
 *  thdu  Thdu *  The generic-table HDU descriptor containing the field.
 *  icol   int    The 1-relative (1->tfields) column number to be
 *                inquired about.
 * Output:
 *  return int    0 - Not a variable length array.
 *                1 - Is a variable length array field.
 */
int iscolvar(Thdu *thdu, int icol)
{
  Bhdu *bhdu;
  int isvar=0;
  if(thdu->type == F_BINTAB) {
    bhdu = (Bhdu *) thdu;
    if(icol>0 && icol<bhdu->tfields)
      isvar = bhdu->fields[icol-1].isvar;
  };
  return isvar;
}

/*.......................................................................
 * Create a new binary-table extension HDU descriptor.
 *
 * Input:
 *  width     int   The number of 8-bit bytes across one row of the table.
 *  nrow      int   The number of rows in the table.
 *  extname  char * A name to refer to this extension by.
 *                  A copy of the string will be made.
 *  extver    int   The version number to assign to the extension. If
 *                  0 - then when the HDU is added to a FITS file, it
 *                  will be assigned the next highest version number for
 *                  name 'extname' and level 'extlevel'.
 *  extlevel  int   The hierachical level of this extension. Send 1 or
 *                  0 if not relevant.
 *  tfields   int   The number of fields per row.
 *  heapsize long   The size to give the heap for variable length arrays
 *                  measured in FITS 8-bit bytes.
 * Output:
 *  return    int   0 - OK.
 *                  1 - Error.
 */ 
Hdu *new_bintab(int nrow, char *extname, int extver, int extlevel, int tfields,
		long heapsize)
{
  Hdu *hdu;
  Bhdu *bhdu;
  int dims[2];
/*
 * Allocate and default initialize a table HDU.
 */
  hdu = new_Hdu(F_BINTAB);
  if(hdu==NULL)
    return hdu;
/*
 * Get the derived version of hdu.
 */
  bhdu = (Bhdu *) hdu;
/*
 * Initialize the HDU base-class fields.
 */
  dims[0] = 0;   /* Lie for the moment since we don't know the table width */
  dims[1] = nrow;
  hdu = ini_Hdu(hdu, B_CHAR, dims, 2, 0, heapsize, 1, extname, extver,
		extlevel, 0, 0);
/*
 * Assign the table base-class fields.
 */
  if(tfields < 1) {
    fprintf(stderr, "new_bintab: Illegal tfields=%d\n", tfields);
    return del_Hdu(hdu);
  };
  bhdu->tfields = tfields;
/*
 * There can be no heap area until the HDU is closed - see end_bhdu.
 */
  bhdu->theap = 0;
/*
 * Allocate the array of field descriptors.
 */
  bhdu->fields = new_Bfields(bhdu);
  if(bhdu->fields==NULL && bhdu->tfields>0)
    return del_Hdu(hdu);
  return hdu;
}

/*.......................................................................
 * Define the structure of a new table column.
 *
 * Input:
 *  hdu      Hdu *  The descriptor of the new table HDU.
 *  icol     int    The 1-relative field number in the table.
 *  tscal double    The scale factor for the new field.
 *  tzero double    The zero offset for the new field.
 *  tform   char *  The mandatory FORTRAN-90 style format for the field.
 *                  The only valid formats are:
 *                   rL (logical), rX (bit), rI (16-bit int), rJ (32-bit int),
 *                   rA (char), rE (float), rD (double), rB (unsigned byte),
 *                   rC (float complex), rM (double complex) and
 *                   rPt (variable length array). Where r is the number
 *                   of elements in the field and t is one of the above
 *                   upper-case type letters except 'P'.
 *  tnull   long    Value used to denote NULL field value - see NOST
 *                  document. Send NONULL if not required.
 *  ttype   char *  An optional name for this field.
 *  tunit   char *  An optional name for the units of measurement of this
 *                  field.
 *  tdisp   char *  An optional string giving the FORTRAN-90 format to
 *                  use to print a field element - see the NOST document.
 *  tdim    char *  An optional means of apportioning the 1-D array in a
 *                  field into an N-dimensional array eg "(10,20,3)".
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int setbfield(Hdu *hdu, int icol, double tscal, double tzero, char *tform,
	      long tnull, char *ttype, char *tunit, char *tdisp, char *tdim)
{
  Bhdu *bhdu = (Bhdu *) hdu;
  Bfield *field;
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL || hdu->type != F_BINTAB) {
    fprintf(stderr, "setbfield: Inappropriate HDU descriptor received\n");
    return 1;
  };
  if(hdu->complete) {
    fprintf(stderr,
	  "setbfield: Illegal attempt to change established HDU structure\n");
    return 1;
  };
  if(tform==NULL) {
    fprintf(stderr, "setbfield: Missing TFORM argument\n");
    return 1;
  };
  if(icol<1 || icol>bhdu->tfields) {
    fprintf(stderr, "setbfield: Column index out of range\n");
    return 1;
  };
  if(tscal==0.0) {
    fprintf(stderr, "setbfield: Error: tscal=0\n");
    return 1;
  };
/*
 * Get the field descriptor.
 */
  field = &bhdu->fields[icol-1];
/*
 * Decompose the tform keyword.
 */
  if(get_format(tform, field))
    return 1;
/*
 * Increment the recorded width of the table to account for the new field.
 */
  hdu->dims[0] += field->width * (field->isvar ? 8 : typesize(field->type));
/*
 * Assign the rest of the members.
 */
  field->tscal = tscal;
  field->tzero = tzero;
  field->tnull = tnull;
  if(!field->tform && tform)
    field->tform = fitsstr(tform);
  if(!field->ttype && ttype)
    field->ttype = fitsstr(ttype);
  if(!field->tunit && tunit)
    field->tunit = fitsstr(tunit);
  if(!field->tdisp && tdisp)
    field->tdisp = fitsstr(tdisp);
  if(!field->tdim && tdim)
    field->tdim = fitsstr(tdim);
  return 0;
}

/*.......................................................................
 * Allocate bhdu->tfields binary-table descriptors.
 *
 * Input:
 *  bhdu     Bhdu *  The descriptor of the binary table HDU in which to
 *                   allocate the array.
 * Output:
 *  return Bfield *  The allocated array of descriptors.
 */
static Bfield *new_Bfields(Bhdu *bhdu)
{
  Bfield *field;
  int i;
/*
 * Sanity checks.
 */
  if(bhdu==NULL) {
    fprintf(stderr, "new_Bfields: Intercepted NULL Bhdu descriptor\n");
    return NULL;
  };
  if(bhdu->fields) {
    fprintf(stderr, "new_Bfields: bhdu->fields allready allocated\n");
    return NULL;
  };
/*
 * Allocate an array of table field descriptors.
 */
  if(bhdu->tfields>0) {
    bhdu->fields = (Bfield *) malloc(sizeof(Bfield) * bhdu->tfields);
    if(bhdu->fields == NULL) {
      fprintf(stderr,
      "new_Bfields: Insufficient memory for binary-table field descriptors\n");
      return NULL;
    };
    for(i=0; i<bhdu->tfields; i++) { /* Initialize members */
      field = &bhdu->fields[i];
      field->type  = DAT_NON;
      field->tbcol = 1;
      field->tscal = 1.0;
      field->tzero = 0.0;
      field->width = 0;
      field->isvar = 0;
      field->form = ' ';
      field->tnull = NONULL;
      field->tform = NULL;
      field->ttype = NULL;
      field->tunit = NULL;
      field->tdisp = NULL;
      field->tdim  = NULL;
    };
  };
  return bhdu->fields;
}

/*.......................................................................
 * Write the header lines for the derived parts of an binary table
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
static ADDFN(add_bhdu)
{
  Bhdu *bhdu = (Bhdu *) hdu;
  Bfield *field;
  size_t tbcol;      /* Bytewize column in table */  
  int i;
/*
 * Check that all table fields have been initialized.
 */
  field = bhdu->fields;
  for(i=0; i<bhdu->tfields; i++, field++) {
    if(field->type == DAT_NON) {
      fprintf(stderr, "add_bhdu: HDU incomplete: Please use setbfield()"
	      " to describe table field %d\n", i+1);
      return 1;
    };
  };
/*
 * Determine field start columns from the widths of each field.
 */
  tbcol = 1;
  field = bhdu->fields;
  for(i=0; i<bhdu->tfields; i++,field++) {
    field->tbcol = tbcol;
    tbcol += field->width * (field->isvar ? 8 : typesize(field->type));
  };
/*
 * Set the heap offset as 0 for now since we do not know how big the
 * HDU will be yet.
 */
  bhdu->theap = 0;
/*
 * Write the Mandatory PCOUNT, GCOUNT and TFIELDS keywords.
 */
  if(wintkey(fits, hdu, "PCOUNT", 0, hdu->pcount, "Random parameter count") ||
     wintkey(fits, hdu, "GCOUNT", 0, hdu->gcount, "Group count") ||
     wintkey(fits, hdu, "TFIELDS",0, bhdu->tfields, "Number of table fields"))
    return 1;
/*
 * Write extension name and version keywords where given.
 */
  if(w_extkeys(fits, hdu))
    return 1;
/*
 * Write the optional THEAP keyword.
 */
  if(bhdu->pcount>0 && wintkey(fits, hdu, "THEAP", 0, bhdu->theap,
			      "Byte offset of heap area"))
    return 1;
/*
 * Write the field description keywords where given.
 */
  if(bhdu->fields) {
    for(i=0; i<bhdu->tfields; i++) {
      field = &bhdu->fields[i];
      if(field->ttype && wstrkey(fits, hdu, "TTYPE", i+1, field->ttype,
				 "Name of this table field"))
	return 1;
      if(field->tunit && wstrkey(fits, hdu, "TUNIT", i+1, field->tunit,
				 "Unit of measurement of this table field"))
	return 1;
      if(field->tnull!=NONULL && wintkey(fits, hdu, "TNULL", i+1, field->tnull,
					"Value used to indicate a NULL item"))
	return 1;
      if(field->tform && wstrkey(fits, hdu, "TFORM", i+1, field->tform,
				 "Format of table field"))
	return 1;
      if(field->tdisp && wstrkey(fits, hdu, "TDISP", i+1, field->tdisp,
				 "Suggested FORTRAN-90 display format"))
	return 1;
      if(field->tdim && wstrkey(fits, hdu, "TDIM", i+1, field->tdim,
				"Dimensions of this field"))
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
 * Copy the derived parts of a Bhdu descriptor into another Bhdu
 * descriptor. This function should be called only via copy_Hdu().
 *
 * Input:
 *  hdu    Hdu * The original hdu to be copied.
 * Output:
 *  return Hdu * The new copied version of 'hdu', or NULL on error.
 */
static COPYFN(cop_bhdu)
{
  Bhdu *old = (Bhdu *) hdu;
  Hdu *new;   /* The new HDU descriptor */
  Bfield *f;  /* A field descriptor from the original HDU descriptor */
  int icol;   /* Column number in table */
/*
 * Create the new ASCII-table descriptor.
 */
  new = new_bintab(old->dims[1], old->extname, 0, old->extlevel, old->tfields,
		   old->pcount);
  if(new==NULL)
    return new;
/*
 * Copy the field descriptors.
 */
  f = old->fields;
  for(icol=1; icol<=old->tfields; icol++, f++) {
    if(setbfield(new, icol, f->tscal, f->tzero, f->tform, f->tnull,
		 f->ttype, f->tunit, f->tdisp, f->tdim))
      return del_Hdu(new);
  };
  return 0;
}

/*.......................................................................
 * Complete the data section of a binary table HDU.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file in which the HDU resides.
 *  hdu      Hdu *  The descriptor of the HDU to be completed.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static ENDFN(end_bhdu)
{
  int saveline; /* The current header line number */
  int waserr=0; /* Error status */
  Bhdu *bhdu = (Bhdu *) hdu;
/*
 * Now that the final number of rows is known, re-write the NAXIS2
 * header keyword.
 */
  saveline = new_hline(hdu, 4);
  waserr=wintkey(fits, hdu, "NAXIS", 2, hdu->dims[1], "Number of table rows.");
  new_hline(hdu, saveline);
/*
 * Deduce the heap offset.
 */
  bhdu->theap = hdu->grpsize - hdu->pcount;
/*
 * If required, update the THEAP keyword value.
 */
  if(hdu->pcount > 0) {
    Fitkey key;   /* The descriptor of the existing GCOUNT keyword */
/*
 * Locate the GCOUNT header line.
 */
    saveline = new_hline(hdu, 0);
/*
 * Find and read the existing GCOUNT header line.
 */
    waserr = waserr || get_key(fits, hdu, "THEAP", DAT_INT, EOH_SEEK, &key);
/*
 * Over-write it if it needs to be changed.
 */
    if(!waserr && KEYINT(key) != bhdu->theap) {
      new_hline(hdu, what_hline(hdu)-1);
      waserr=wintkey(fits, hdu, "THEAP", 0, bhdu->theap, "Number of groups.");
    };
/*
 * Restore the header-line pointer.
 */
    new_hline(hdu, saveline);
  };
  return waserr!=0;
}
