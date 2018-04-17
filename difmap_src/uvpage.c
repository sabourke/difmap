#include <stdlib.h>
#include <stdio.h>
 
#include "logio.h"
#include "recio.h"
#include "uvpage.h"

/*
 * Define a Mvis structure used to initialize visibilities in the
 * buffer.
 */
static Mvis zero_vis = {0.0f,0.0f};

static UVpage *uvpmemerr(UVpage *uvp);

/*.......................................................................
 * Open a uvmodel.scr scratch file and return a descriptor to be used
 * in I/O to the file.
 *
 * Input:
 *  ntime     int    The number of integrations in each model.
 *  nbase     int    The number of baselines in each model.
 *  nif       int    The number of IFs for which a UV model is to be stored.
 * Output:
 *  return UVpage *  The initialized descriptor, or NULL on error.
 */
UVpage *new_UVpage(int ntime, int nbase, int nif)
{
  UVpage *uvp;   /* The descriptor to be returned */
/*
 * Check the validity of the arguments.
 */
  if(ntime <= 0 || nbase <= 0 || nif <= 0) {
    lprintf(stderr,
	    "new_UVpage: Arguments specify a -ve number of visibilities.\n");
    return NULL;
  };
/*
 * Allocate the container.
 */
  uvp = malloc(sizeof(UVpage));
  if(uvp==NULL)
    return uvpmemerr(uvp);
/*
 * Initialize the members of the descriptor at least to the extent that
 * the descriptor can subsequently be passed to del_UVpage().
 */
  uvp->rio = NULL;
  uvp->ntime = ntime;
  uvp->nbase = nbase;
  uvp->nif = nif;
  uvp->mvis = NULL;
  uvp->ioerr = 0;
/*
 * Open the binary scratch file, using one IF as the record length.
 */
  uvp->rio = new_Recio("uvmodel.scr", IS_SCR, 0, ntime * nbase * sizeof(Mvis));
  if(uvp->rio == NULL)
    return del_UVpage(uvp);
/*
 * Allocate a buffer to be used when reading and writing to the
 * uvmodel.scr.
 */
  uvp->mvis = malloc(nbase * sizeof(Mvis));
  if(uvp->mvis == NULL)
    return uvpmemerr(uvp);
/*
 * Return the initialized container.
 */
  return uvp;
}

/*.......................................................................
 * Private cleanup function of new_UVpage() for memory allocation failures.
 *
 * Input:
 *  uvp     UVpage *   The partially initialized UVpage descriptor.
 * Output:
 *  return  UVpage *   Allways NULL.
 */
static UVpage *uvpmemerr(UVpage *uvp)
{
  lprintf(stderr, "new_UVpage: Insufficient memory.\n");
  return del_UVpage(uvp);
}

/*.......................................................................
 * Close and delete a uvmodel.scr paging file.
 *
 * Input:
 *  uvp     UVpage * The descriptor to be deleted.
 * Output:
 *  return  UVpage * Allways NULL. Use like uvp=del_UVpage(uvp);
 */
UVpage *del_UVpage(UVpage *uvp)
{
  if(uvp) {
/*
 * Close and delete the uvmodel.scr file.
 */
    if(uvp->rio)
      del_Recio(uvp->rio);
/*
 * Delete the visbility I/O buffer.
 */
    if(uvp->mvis)
      free(uvp->mvis);
  };
  return NULL;
}

/*.......................................................................
 * Read an integration worth of the model visibilities of a given IF.
 * The visibilities will be stored in uvp->mvis[0..uvp->nbase-1].
 *
 * Input:
 *  uvp     UVpage *   The paging descriptor.
 *  ut         int     The index of the integration to be read.
 *  cif        int     The index of the IF to be read from.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
int uvp_read(UVpage *uvp, int ut, int cif)
{
/*
 * Check validity of uvp.
 */
  if(uvp_error(uvp, "uvp_read"))
    return 1;
/*
 * Check integration index.
 */
  if(ut < 0 || ut >= uvp->ntime) {
    lprintf(stderr, "uvp_read: Integration index out of range.\n");
    return 1;
  };
/*
 * IF integration index.
 */
  if(cif<0 || cif>=uvp->nif)  {
    lprintf(stderr, "uvp_read: IF index out of range.\n");
    return 1;
  };
/*
 * Position the file if necessary.
 */
  if(rec_seek(uvp->rio, cif, ut * uvp->nbase * sizeof(Mvis))) {
    uvp->ioerr = 1;
    return 1;
  };
/*
 * Read from the scratch file into the buffer.
 */
  if(rec_read(uvp->rio, uvp->nbase, sizeof(Mvis), uvp->mvis) < uvp->nbase) {
    lprintf(stderr, "uvp_read: Error reading from scratch file.\n");
    uvp->ioerr = 1;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Write an integration worth of the model visibilities of a given IF.
 * The visibilities will be extracted from uvp->mvis[0..uvp->nbase-1].
 *
 * Input:
 *  uvp     UVpage *   The paging descriptor.
 *  ut         int     The index of the integration to be written to.
 *  cif        int     The index of the IF to be written in.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
int uvp_write(UVpage *uvp, int ut, int cif)
{
/*
 * Check validity of uvp.
 */
  if(uvp_error(uvp, "uvp_write"))
    return 1;
/*
 * Check integration index.
 */
  if(ut < 0 || ut >= uvp->ntime) {
    lprintf(stderr, "uvp_write: Integration index out of range.\n");
    return 1;
  };
/*
 * IF integration index.
 */
  if(cif<0 || cif>=uvp->nif)  {
    lprintf(stderr, "uvp_write: IF index out of range.\n");
    return 1;
  };
/*
 * Position the file if necessary.
 */
  if(rec_seek(uvp->rio, cif, ut * uvp->nbase * sizeof(Mvis))) {
    uvp->ioerr = 1;
    return 1;
  };
/*
 * Read from the scratch file into the buffer.
 */
  if(rec_write(uvp->rio, uvp->nbase, sizeof(Mvis), uvp->mvis) < uvp->nbase) {
    lprintf(stderr, "uvp_write: Error reading from scratch file.\n");
    uvp->ioerr = 1;
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Check the validity of the UVpage descriptor.
 *
 * Input:
 *  uvp   UVpage *   The UV model paging file descriptor.
 *  fname   char *   The name of the calling function.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int uvp_error(UVpage *uvp, const char *fname)
{
/*
 * No descriptor?
 */
  if(uvp==NULL) {
    lprintf(stderr, "%s: Intercepted NULL UVpage descriptor.\n", fname);
    return 1;
  };
/*
 * Previous I/O error?
 */
  if(uvp->ioerr)
    return 1;
  return 0;
}

/*.......................................................................
 * Clear the whole output buffer.
 *
 * Input:
 *  uvp  UVpage *  The UV model paging file descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int uvp_clear(UVpage *uvp)
{
  Mvis *mvis;   /* Pointer to array of uvp->nbase visibilities */
  long nvis;    /* Number of visibilities left to be cleared */
/*
 * Check the validity of uvp.
 */
  if(uvp_error(uvp, "uvp_clear"))
    return 1;
/*
 * Clear the buffer.
 */
  nvis = uvp->nbase;
  mvis = uvp->mvis;
  while(nvis--)
    *mvis++ = zero_vis;
  return 0;
}

/*.......................................................................
 * Make sure that a UV paging file is up to date by flushing all
 * I/O.
 *
 * Input:
 *  uvp  UVpage *  The UV paging file descriptor.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int uvp_flush(UVpage *uvp)
{
  return uvp ? rec_flush(uvp->rio) : 0;
}
