#ifndef ifpage_h
#define ifpage_h

#include <stdlib.h>

/* Declare a struct to hold one visibility in amp,phase form. */

typedef struct {
  float amp;      /* Amplitude of visibility */
  float phs;      /* Phase of visibility */
  float wt;       /* Visibility weight (1/variance) */
} Dvis;

/* 
 * Declare the descriptor used to read and write data to and from the
 * data paging file.
 */
typedef struct {
  Recio *rio;    /* Record I/O descriptor */
  int ioerr;     /* True after a record I/O error */
  Dvis *dvis;    /* A buffer of sufficient size to contain one integration */
  long first;    /* The index of the first visibility in the buffer */
  long nread;    /* The number of visibilities in the buffer */
  int nbase;     /* The number of baseline visibilities in dvis[] */
  int ntime;     /* The number of integrations in the file */
  int nif;       /* The number  of IF's in the file */
  int cif;       /* Index of currently selected IF */
} IFpage;

/* Open a new binary ifdata.scr scratch file */

IFpage *new_IFpage(int nif, int nbase, int ntime);

/* Close and delete an ifdata.scr scratch file */

IFpage *del_IFpage(IFpage *ip);

/* Read data of integration ut from range set by ip_range() */

int ip_read(IFpage *ip, long ut);

/* Write data to integration ut from range set by ip_range() */

int ip_write(IFpage *ip, long ut);

/* Check validity of an IFpage descriptor */

int ip_error(IFpage *ip, const char *fname);

/*----------------------------------------------------------------------
 * The following functions may legally be called with ip==NULL, in which
 * case they quietly do nothing and return as though successfull.
 */

/* Set subsequent integration read-write range */

int ip_range(IFpage *ip, int ifa, int ba, int bb);

/* Clear the whole IFpage integration buffer */

int ip_clear(IFpage *ip);

/* Flush pending I/O to the paging file */

int ip_flush(IFpage *ip);

#endif
