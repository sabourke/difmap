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
#include "thdu.h"
#include "fitkey.h"

/* Virtual function table for each known table type */

extern Tabfn atabfn;
extern Tabfn btabfn;

static struct Tabtab {
  Hdutype type;
  Tabfn *fns;
} tabtab[]={
  {F_TABLE,  &atabfn},  /* ASCII tables */
  {F_BINTAB, &btabfn},  /* Binary tables */
};

static Tabfn *whatthdu(Thdu *thdu);

/*.......................................................................
 * Return the virtual function table entry for a given table type.
 *
 * Input:
 *  thdu     Thdu * The generic table HDU descriptor. If this turns out
 *                  not to be a table HDU or if thdu==NULL then an error
 *                  message will be issued and NULL returned.
 * Output:
 *  return  Tabfn * The table entry - or NULL on error.
 */
static Tabfn *whatthdu(Thdu *thdu)
{
  int i;
/*
 * Valid descriptor?
 */
  if(thdu==NULL) {
    fprintf(stderr, "whatthdu: NULL Thdu descriptor intercepted\n");
    return NULL;
  };
/*
 * Search for the HDU type in the subset contained in tabtab.
 */
  for(i=0; i<sizeof(tabtab)/sizeof(struct Tabtab); i++) {
    if(thdu->type == tabtab[i].type)
      return tabtab[i].fns;
  };
  fprintf(stderr, "whatthdu: HDU does not point to a table HDU\n");
  return NULL;
}

/*.......................................................................
 * Clear the derived part of an generic-table header descriptor to the
 * point at which it can be safely sent to del_Hdu(). This should be
 * called only by the equivalent function for the specific table type.
 *
 * Input:
 *  hdu     Hdu *  The base-class pointer to a Thdu descriptor.
 * Output:
 *  return  Hdu *  'hdu', or NULL if an error occurs and hdu is deleted.
 */
NEWFN(new_table)
{
  Thdu *thdu = (Thdu *) hdu;
  thdu->tfields = 0;
  return hdu;
}

/*.......................................................................
 * Delete the derived parts of a generic table HDU descriptor.
 *
 * Input:
 *  hdu      Hdu *  The base-class pointer to the Thdu descriptor.
 */
DELFN(del_table)
{
/*
 * There is nothing to delete at the moment - this function is here
 * solely for the purposes of future expansion.
 */
  return;
}

/*.......................................................................
 * Find the HDU of a table by its name and version number.
 *
 * Input:
 *  fits     Fits *  The FITS descriptor.
 *  extname  char *  The name of the extension to be found. Trailing
 *                   spaces are ignored in the comparisons.
 *  extver    int    The version number of the table, or 0 for the last
 *                   table of this type in the FITS file.
 *  prev      Hdu *  The HDU to search forwards from. Send NULL to
 *                   start from the start HDU in the file. The first
 *                   hdu checked is prev->next.
 * Output:
 *  return   Thdu *  The pointer to the generic table HDU descriptor
 *                   or NULL if not found.
 */
Thdu *find_table(Fits *fits, char *extname, int extver, Hdu *prev)
{
  return (Thdu *) find_hdu(fits, F_TABLE | F_BINTAB, extname, extver, prev);
}

/*.......................................................................
 * Return the number of the table-column that has a specified name.
 *
 * Input:
 *  thdu    Thdu *  The generic table HDU descriptor.
 *  ttype   char *  The name of the column to be sought - Only as many
 *                  characters as are given will be compared.
 *  fixlen   int    If > 0, then this defines the max number of characters
 *                  to be compared.
 * Output:
 *  return   int    1-relative column number, or 0 if not found.
 */
COL_FINDFN(find_column)
{
  Tabfn *fntab;
/*
 * Delegate the task to a method function of the appropriate
 * table type.
 */
  fntab = whatthdu(thdu);
  return fntab ? (*fntab->findfn)(thdu, ttype, fixlen) : 0;
}

/*.......................................................................
 * Return the data-type of a given column.
 *
 * Input:
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 * Output:
 *  return  Fittype    The data-type of the column. On error DAT_NON==0 is
 *                     returned.
 */
COL_TYPEFN(col_type)
{
  Tabfn *fntab;
/*
 * Sanity check the HDU type.
 */
  if(thdu==NULL) {
    fprintf(stderr, "col_type: NULL HDU descriptor received\n");
    return DAT_NON;
  };
  if(icol < 1 || icol > thdu->tfields) {
    fprintf(stderr, "col_type: Out-of-range column index received\n");
    return DAT_NON;
  };
/*
 * Delegate the task to a method function of the appropriate
 * table type.
 */
  fntab = whatthdu(thdu);
  return fntab ? (*fntab->typefn)(thdu, icol) : DAT_NON;
}

/*.......................................................................
 * Return the dimension of entries at a given row of a given column.
 *
 * Input:
 *  fits       Fits *  The FITS file descriptor.
 *  thdu       Thdu *  The table HDU-descriptor.
 *  icol        int    The 1-relative number of the column.
 *  irow        int    The 1-relative row number. This is ignored except
 *                     for columns containing variable length arrays.
 * Output:
 *  return      int    The number of elements in the specified field.
 *                     This may be 0. 0 is also returned on error.
 */
COL_DIMFN(col_dim)
{
  Tabfn *fntab;
/*
 * Sanity check the HDU type.
 */
  if(thdu==NULL) {
    fprintf(stderr, "col_dim: NULL HDU descriptor received\n");
    return DAT_NON;
  };
  if(icol < 1 || icol > thdu->tfields) {
    fprintf(stderr, "col_dim: Out-of-range column index received\n");
    return DAT_NON;
  };
/*
 * Delegate the task to a method function of the appropriate
 * table type.
 */
  fntab = whatthdu(thdu);
  return fntab ? (*fntab->dimfn)(fits, thdu, icol, irow) : 0;
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
 *                     conversions will be performed between the type
 *                     stored in the FITS file and this type - where this
 *                     is meaningful.
 *  doscale     int    If true apply offset and scale factors.
 *  flags  Fitsflag *  If not NULL, then this should point to an array
 *                     of nobj elements to be used to record which elements
 *                     of data[] are blanked 1=blanked, 0=OK.
 *  first       int    The 0-relative index of the first element to get.
 *  ndata       int    The max number of elements to return from the field.
 *                     NB. String values will not be '\0' terminated.
 * Output:
 *  data       void *  Supply an array of the type encoded in 'type' and
 *                     having at least 'ndata' elements.
 *  return     long    Number of objects read. This will be less than
 *                     ndata if the field contains fewer than ndata objects.
 *                     On error 0 is returned.
 */
COL_VALFN(rcolumn)
{
  Tabfn *fntab;
/*
 * Sanity check the HDU type.
 */
  if(thdu==NULL) {
    fprintf(stderr, "rcolumn: NULL HDU descriptor received\n");
    return 0;
  };
  if(icol < 1 || icol > thdu->tfields) {
    fprintf(stderr, "rcolumn: Out-of-range column index received\n");
    return 0;
  };
/*
 * The row index may exceed the currently recorded number of rows only
 * if the data section is marked as incomplete.
 */
  if(irow < 1 || (thdu->state!=HDU_DATA && irow > thdu->dims[1])) {
    fprintf(stderr, "rcolumn: Out of range (1-%d) row index (%d) rejected.\n",
	    thdu->dims[1], irow);
    return 0;
  };
/*
 * Delegate the task to a method function of the appropriate
 * table type.
 */
  fntab = whatthdu(thdu);
  return fntab ? (*fntab->valfn)(fits, thdu, icol, irow, type, doscale, flags, first, ndata, data) : 0;
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
 *                     ndata if the field contains fewer than ndata objects.
 *                     On error 0 is returned.
 */
COL_SETFN(wcolumn)
{
  Tabfn *fntab;
/*
 * Sanity check the HDU type.
 */
  if(thdu==NULL) {
    fprintf(stderr, "wcolumn: NULL HDU descriptor received\n");
    return 0;
  };
  if(icol < 1 || icol > thdu->tfields) {
    fprintf(stderr, "wcolumn: Out-of-range column index received\n");
    return 0;
  };
  if(irow < 1) {
    fprintf(stderr, "wcolumn: Illegal negative or 0 row index rejected.\n");
    return 0;
  };
/*
 * The row index may exceed the currently recorded number of rows only
 * if the data section is marked as incomplete.
 */
  if(irow > thdu->dims[1]) {
   if(thdu->state == HDU_DATA) {
     thdu->dims[1] = irow;      /* Expand the recorded number of rows */
   } else {
     fprintf(stderr, "wcolumn: Can't expand table to accomodate row %d.\n",
	     irow);
     return 0;
   };
  };
/*
 * Delegate the task to a method function of the appropriate
 * table type.
 */
  fntab = whatthdu(thdu);
  return fntab ? (*fntab->setfn)(fits, thdu, icol, irow, type, doscale, flags, first, ndata, data) : 0;
}

/*.......................................................................
 * Return the name of a given column.
 *
 * Input:
 *  thdu     Thdu *  The generic table descriptor.
 *  icol      int    The number of the column to be looked up.
 * Output:
 *  return   char *  The column name or NULL if not known or on error.
 */
COL_NAMEFN(col_name)
{
  Tabfn *fntab;   /* Vitual function table entry for table */
/*
 * Sanity check the HDU type.
 */
  if(thdu==NULL) {
    fprintf(stderr, "col_name: NULL HDU descriptor received\n");
    return NULL;
  };
  if(icol < 1 || icol > thdu->tfields) {
    fprintf(stderr, "col_name: Column index (%d) out-of-range\n", icol);
    return NULL;
  };
/*
 * Delegate the task to a method function of the appropriate
 * table type.
 */
  fntab = whatthdu(thdu);
  return fntab ? (*fntab->namefn)(thdu, icol) : NULL;
}

/*.......................................................................
 * Public function to return the number of rows in a table.
 *
 * Input:
 *  thdu   Thdu *  The generic table HDU descriptor of the table.
 * Output:
 *  return  int    The number of rows in the table, or 0 on error.
 */
int numrow(Thdu *thdu)
{
  if(thdu->type != F_BINTAB && thdu->type != F_TABLE) {
    fprintf(stderr, "numrow: The HDU descriptor is not of a table\n");
    return 0;
  };
  return thdu->dims[1];
}

/*.......................................................................
 * Public function to return the number of fields in a table.
 *
 * Input:
 *  thdu   Thdu *  The generic table HDU descriptor of the table.
 * Output:
 *  return  int    The number of fields in the table, or 0 on error.
 */
int numcol(Thdu *thdu)
{
  if(thdu->type != F_BINTAB && thdu->type != F_TABLE) {
    fprintf(stderr, "numcol: The HDU descriptor is not of a table\n");
    return 0;
  };
  return thdu->tfields;
}
