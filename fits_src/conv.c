#include <stdio.h>
#include "fits.h"

/* Declare the form of a type conversion function */

#define CONV_FN(fn) void (fn)(long ndata, void *adata, double zero, double scal, void *bdata)
typedef CONV_FN(Conv_fn);

/* List type-conversion functions */

static Conv_fn b2b;
static Conv_fn b2sh;
static Conv_fn b2i;
static Conv_fn b2lg;
static Conv_fn b2f;
static Conv_fn b2d;
static Conv_fn sh2b;
static Conv_fn sh2sh;
static Conv_fn sh2i;
static Conv_fn sh2lg;
static Conv_fn sh2f;
static Conv_fn sh2d;
static Conv_fn i2b;
static Conv_fn i2sh;
static Conv_fn i2i;
static Conv_fn i2lg;
static Conv_fn i2f;
static Conv_fn i2d;
static Conv_fn i2l;
static Conv_fn lg2b;
static Conv_fn lg2sh;
static Conv_fn lg2i;
static Conv_fn lg2lg;
static Conv_fn lg2f;
static Conv_fn lg2d;
static Conv_fn f2b;
static Conv_fn f2sh;
static Conv_fn f2i;
static Conv_fn f2lg;
static Conv_fn f2f;
static Conv_fn f2d;
static Conv_fn d2b;
static Conv_fn d2sh;
static Conv_fn d2i;
static Conv_fn d2lg;
static Conv_fn d2f;
static Conv_fn d2d;
static Conv_fn c2c;
static Conv_fn l2l;
static Conv_fn l2i;
static Conv_fn sc2sc;
static Conv_fn sc2dc;
static Conv_fn dc2sc;
static Conv_fn dc2dc;
static Conv_fn st2st;

/* Create a table of conversion functions. Where a conversion is not handled */

static Conv_fn *conv_table[13][13] = {
/*           SHT, INT,  LNG, FLT, DBL,CHR,BYT,BIT,LOG, SCMP, DCMP,  COM,  STR */
/* SHT */ {sh2sh,sh2i,sh2lg,sh2f,sh2d, 0,sh2b,  0,  0,    0,    0,    0,    0},
/* INT */ { i2sh, i2i, i2lg, i2f, i2d, 0, i2b,  0,i2l,    0,    0,    0,    0},
/* LNG */ {lg2sh,lg2i,lg2lg,lg2f,lg2d, 0,lg2b,  0,  0,    0,    0,    0,    0},
/* FLT */ { f2sh, f2i, f2lg, f2f, f2d, 0, f2b,  0,  0,    0,    0,    0,    0},
/* DBL */ { d2sh, d2i, d2lg, d2f, d2d, 0, d2b,  0,  0,    0,    0,    0,    0},
/* CHR */ {    0,   0,    0,   0,   0,c2c,  0,  0,  0,    0,    0,    0,    0},
/* BYT */ { b2sh, b2i, b2lg, b2f, b2d, 0, b2b,  0,  0,    0,    0,    0,    0},
/* BIT */ {    0,   0,    0,   0,   0, 0,   0,l2l,  0,    0,    0,    0,    0},
/* LOG */ {    0, l2i,    0,   0,   0, 0,   0,  0,l2l,    0,    0,    0,    0},
/* SCMP */{    0,   0,    0,   0,   0, 0,   0,  0,  0,sc2sc,sc2dc,    0,    0},
/* DCMP */{    0,   0,    0,   0,   0, 0,   0,  0,  0,dc2sc,dc2dc,    0,    0},
/* COM */ {    0,   0,    0,   0,   0, 0,   0,  0,  0,    0,    0,st2st,st2st},
/* STR */ {    0,   0,    0,   0,   0, 0,   0,  0,  0,    0,    0,st2st,st2st},
};

/*.......................................................................
 * Convert an array of one data-type to one of another.
 *
 * Input:
 *  ndata      long    The number of elements in 'adata[]' and 'bdata[]'.
 *  atype   Fittype    The data-type of the input array.
 *  adata      void *  The input array of dimension 'ndata' and type 'atype'.
 *  zero     double    Arithmetic data-types are offset by this value
 *                     during conversion.
 *  scal     double    Arithmetic data-types are scaled by this value
 *                     during conversion.
 *  btype   Fittype    The data-type of the output array (bdata).
 * Output:
 *  bdata      void *  The output array of dimension 'ndata' and type 'btype'.
 *  return      int    0 - OK.
 *                     1 - Error, eg. Unhandled conversion.
 */
int typeconv(long ndata, Fittype atype, void *adata, double zero, double scal,
	     Fittype btype, void *bdata)
{
  Conv_fn *convfn=0;  /* Pointer to conversion function */
/*
 * Work out the indexes into the function table for the given types.
 */
  int aslot = ((int)atype) - 1;
  int bslot = ((int)btype) - 1;
/*
 * Check that the types are in range of the lookup table.
 */
  if(aslot>=0 && aslot<sizeof(conv_table[0]) &&
     bslot>=0 && bslot<sizeof(conv_table[0]))
     convfn = conv_table[aslot][bslot];
/*
 * Conversion function available?
 */
  if(convfn)
    convfn(ndata, adata, zero, scal, bdata);
  else {
    fprintf(stderr, "typeconv: Unhandled conversion from (%s) to (%s)\n",
	    typename(atype), typename(btype));
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Convert a single value of one data-type to one of another.
 *
 * Input:
 *  atype   Fittype    The data-type of the input array.
 *  adata      void *  The input array of dimension 'ndata' and type 'atype'.
 *  zero     double    Arithmetic data-types are offset by this value
 *                     during conversion.
 *  scal     double    Arithmetic data-types are scaled by this value
 *                     during conversion.
 *  btype   Fittype    The data-type of the output array (bdata).
 * Output:
 *  bdata      void *  The output array of dimension 'ndata' and type 'btype'.
 *  return      int    0 - OK.
 *                     1 - Error, eg. Unhandled conversion.
 */
int stypeconv(Fittype atype, void *adata, double zero, double scal,
	      Fittype btype, void *bdata)
{
  Conv_fn *convfn=0;  /* Pointer to conversion function */
/*
 * Work out the indexes into the function table for the given types.
 */
  int aslot = ((int)atype) - 1;
  int bslot = ((int)btype) - 1;
/*
 * Check that the types are in range of the lookup table.
 */
  if(aslot>=0 && aslot<sizeof(conv_table[0]) &&
     bslot>=0 && bslot<sizeof(conv_table[0]))
     convfn = conv_table[aslot][bslot];
/*
 * Conversion function available?
 */
  if(convfn)
    convfn(1L, adata, zero, scal, bdata);
  else {
    fprintf(stderr, "typeconv: Unhandled conversion from (%s) to (%s)\n",
	    typename(atype), typename(btype));
    return 1;
  };
  return 0;
}

/*-----------------------------------------------------------------------
 * All the following functions copy a scaled and offset version of a
 * specific type of input array to a specific type of output array.
 *
 * All take the folowing arguments.
 *
 * Input:
 *  ndata      long    The number of elements in 'adata[]' and 'bdata[]'.
 *  adata      void *  A pointer to the input array array of dimension
 *                     'ndata' after being cast to (void *).
 *  zero     double    Arithmetic data-types are offset by this value
 *                     during conversion.
 *  scal     double    Arithmetic data-types are scaled by this value
 *                     during conversion.
 * Output:
 *  bdata      void *  A pointer to the output array array of dimension
 *                     'ndata' after being cast to (void *).
 *-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------
 *  DAT_BYT to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_BYT -> DAT_BYT */

static CONV_FN(b2b)
{
  long i;
  unsigned char *aptr = adata;
  unsigned char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_BYT -> DAT_SHT */

static CONV_FN(b2sh)
{
  long i;
  unsigned char *aptr = adata;
  short *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_BYT -> DAT_INT */

static CONV_FN(b2i)
{
  long i;
  unsigned char *aptr = adata;
  int *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_BYT -> DAT_LNG */

static CONV_FN(b2lg)
{
  long i;
  unsigned char *aptr = adata;
  long *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_BYT -> DAT_FLT */

static CONV_FN(b2f)
{
  long i;
  unsigned char *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_BYT -> DAT_DBL */

static CONV_FN(b2d)
{
  long i;
  unsigned char *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_SHT to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_SHT -> DAT_BYT */

static CONV_FN(sh2b)
{
  long i;
  short *aptr = adata;
  unsigned char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_SHT -> DAT_SHT */

static CONV_FN(sh2sh)
{
  long i;
  short *aptr = adata;
  short *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_SHT -> DAT_INT */

static CONV_FN(sh2i)
{
  long i;
  short *aptr = adata;
  int *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_SHT -> DAT_LNG */

static CONV_FN(sh2lg)
{
  long i;
  short *aptr = adata;
  long *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_SHT -> DAT_FLT */

static CONV_FN(sh2f)
{
  long i;
  short *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_SHT -> DAT_DBL */

static CONV_FN(sh2d)
{
  long i;
  short *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_INT to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_INT -> DAT_BYT */

static CONV_FN(i2b)
{
  long i;
  int *aptr = adata;
  unsigned char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_INT -> DAT_SHT */

static CONV_FN(i2sh)
{
  long i;
  int *aptr = adata;
  short *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_INT -> DAT_INT */

static CONV_FN(i2i)
{
  long i;
  int *aptr = adata;
  int *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_INT -> DAT_LNG */

static CONV_FN(i2lg)
{
  long i;
  int *aptr = adata;
  long *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_INT -> DAT_FLT */

static CONV_FN(i2f)
{
  long i;
  int *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_INT -> DAT_DBL */

static CONV_FN(i2d)
{
  long i;
  int *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_INT -> DAT_LOG */

static CONV_FN(i2l)
{
  long i;
  int *aptr = adata;
  char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = aptr[i] ? 'T' : 'F';
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_LNG to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_LNG -> DAT_BYT */

static CONV_FN(lg2b)
{
  long i;
  long *aptr = adata;
  unsigned char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_LNG -> DAT_SHT */

static CONV_FN(lg2sh)
{
  long i;
  long *aptr = adata;
  short *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_LNG -> DAT_INT */

static CONV_FN(lg2i)
{
  long i;
  long *aptr = adata;
  int *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_LNG -> DAT_LNG */

static CONV_FN(lg2lg)
{
  long i;
  long *aptr = adata;
  long *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_LNG -> DAT_FLT */

static CONV_FN(lg2f)
{
  long i;
  long *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_LNG -> DAT_DBL */

static CONV_FN(lg2d)
{
  long i;
  long *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_FLT to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_FLT -> DAT_BYT */

static CONV_FN(f2b)
{
  long i;
  float *aptr = adata;
  unsigned char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_FLT -> DAT_SHT */

static CONV_FN(f2sh)
{
  long i;
  float *aptr = adata;
  short *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_FLT -> DAT_INT */

static CONV_FN(f2i)
{
  long i;
  float *aptr = adata;
  int *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_FLT -> DAT_LNG */

static CONV_FN(f2lg)
{
  long i;
  float *aptr = adata;
  long *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_FLT -> DAT_FLT */

static CONV_FN(f2f)
{
  long i;
  float *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_FLT -> DAT_DBL */

static CONV_FN(f2d)
{
  long i;
  float *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_DBL to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_DBL -> DAT_BYT */

static CONV_FN(d2b)
{
  long i;
  double *aptr = adata;
  unsigned char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_DBL -> DAT_SHT */

static CONV_FN(d2sh)
{
  long i;
  double *aptr = adata;
  short *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_DBL -> DAT_INT */

static CONV_FN(d2i)
{
  long i;
  double *aptr = adata;
  int *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_DBL -> DAT_LNG */

static CONV_FN(d2lg)
{
  long i;
  double *aptr = adata;
  long *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_DBL -> DAT_FLT */

static CONV_FN(d2f)
{
  long i;
  double *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/* DAT_DBL -> DAT_DBL */

static CONV_FN(d2d)
{
  long i;
  double *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_CHR to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_CHR -> DAT_CHR */

static CONV_FN(c2c)
{
  long i;
  char *aptr = adata;
  char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = zero + scal * aptr[i];
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_BIT to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_BIT -> DAT_BIT  [use l2l()] */


/*-----------------------------------------------------------------------
 *  DAT_LOG to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_LOG -> DAT_LOG */

static CONV_FN(l2l)
{
  long i;
  char *aptr = adata;
  char *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = aptr[i];
  return;
}

/* DAT_LOG -> DAT_INT */

static CONV_FN(l2i)
{
  long i;
  char *aptr = adata;
  int *bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = aptr[i]=='T';
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_SCMP to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_SCMP -> DAT_SCMP */

static CONV_FN(sc2sc)
{
  long i;
  float *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i+=2) {
    bptr[i] = zero + scal * aptr[i];
    bptr[i+1] = scal * aptr[i+1];
  };
  return;
}

/* DAT_SCMP -> DAT_DCMP */

static CONV_FN(sc2dc)
{
  long i;
  float *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i+=2) {
    bptr[i] = zero + scal * aptr[i];
    bptr[i+1] = scal * aptr[i+1];
  };
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_DCMP to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_DCMP -> DAT_SCMP */

static CONV_FN(dc2sc)
{
  long i;
  double *aptr = adata;
  float *bptr = bdata;
  for(i=0; i<ndata; i+=2) {
    bptr[i] = zero + scal * aptr[i];
    bptr[i+1] = scal * aptr[i+1];
  };
  return;
}

/* DAT_DCMP -> DAT_DCMP */

static CONV_FN(dc2dc)
{
  long i;
  double *aptr = adata;
  double *bptr = bdata;
  for(i=0; i<ndata; i+=2) {
    bptr[i] = zero + scal * aptr[i];
    bptr[i+1] = scal * aptr[i+1];
  };
  return;
}

/*-----------------------------------------------------------------------
 *  DAT_COM to whatever.
 *-----------------------------------------------------------------------*/

/* DAT_STR || DAT_COM -> DAT_STR || DAT_COM */

static CONV_FN(st2st)
{
  long i;
  char **aptr = adata;
  char **bptr = bdata;
  for(i=0; i<ndata; i++)
    bptr[i] = aptr[i];
  return;
}
