#ifndef dpage_h
#define dpage_h

#include <stdlib.h>

/* Declare a struct to hold one complex visibility as real,imaginary,weight */

typedef struct {
  float re;        /* Real part of visibility */
  float im;        /* Imaginary part of visibility */
  float wt;        /* Weight of visibility (1/variance) */
} Cvis;

/*
 * Declare a hierarchical tree to refer to the visibilities.
 * This tree currently reflects the order of visibilities in the
 * integration record.
 */

typedef struct {  /* Each baseline cites dp->npol visibilities */
  Cvis *pol;
} Dbase;

typedef struct {  /* Each spectral-line channel cites dp->nbase baselines */ 
  Dbase *base;
} Dchan;

typedef struct {  /* Each IF cites dp->nchan spectral-line channels */
  Dchan *chan;
} Dif;

/* 
 * Declare the descriptor used to read and write data to and from the
 * data paging file.
 */
typedef struct {
  Recio *rio;    /* Record I/O descriptor */
  int ioerr;     /* True after a record I/O error */
  int ut;        /* The UT currently in cvis[], or -1 if not initialized */
  Cvis *cvis;    /* A buffer of sufficient size to contain one integration */
  size_t nvis;   /* Number of visibilities per integration */
  Dif *ifs;      /* Pointer to array of nif IFs containers. */
  long first;    /* The index of the first visibility in the buffer */
  long nbuff;    /* The number of visibilities in the buffer */
  int nbase;     /* The number of baselines in the file */
  int ntime;     /* The number of integrations in the file */
  int nchan;     /* The number of spectral-line channels in the file */
  int nif;       /* The number  of IF's in the file */
  int npol;      /* The number of stokes parameters */
  long soff;     /* Indexing offset between stokes in cvis[] */
  long boff;     /* Indexing offset between baselines in cvis[] */
  long coff;     /* Indexing offset between frequencies in cvis[] */
  long ioff;     /* Indexing offset between IFs in cvis[] */
  int ca,cb;     /* Indexes of first and last channels to be transferred */
  int ia,ib;     /* Indexes of first and last IFs to be transferred */
  int sa,sb;     /* Indexes of first and last stokes to be transferred */
  int ba,bb;     /* Indexes of first and last baselines to be transferred */
} Dpage;

/* Open a uvdata.scr paging file */

Dpage *new_Dpage(int ntime, int nbase, int nchan, int nif, int npol);

/* Close and delete a uvdata.scr paging file */

Dpage *del_Dpage(Dpage *dp);

/* Read data of integration ut from range set by dp_range() */

int dp_read(Dpage *dp, long ut);

/* Write data to integration ut from range set by dp_range() */

int dp_write(Dpage *dp, long ut);

/*----------------------------------------------------------------------
 * The following functions may legally be called with dp==NULL, in which
 * case they quietly do nothing and return as though successfull.
 */

/* Set subsequent integration read-write range */

int dp_crange(Dpage *dp, int ca, int cb);
int dp_irange(Dpage *dp, int ifa, int ifb);
int dp_brange(Dpage *dp, int ba, int bb);
int dp_srange(Dpage *dp, int sa, int sb);

/* Clear the whole dpage integration buffer for a new integration */

int dp_clear(Dpage *dp, int ut);

/* Flush pending I/O to the paging file */

int dp_flush(Dpage *dp);

#endif
