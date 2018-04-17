#include <stdio.h>
#include <stdlib.h>

#include "logio.h"
#include "obs.h"

/*-----------------------------------------------------------------------
 * Handle allocation, re-sizing, default initialization, and deletion of
 * Binan AIPS binary AN table descriptors.
 *----------------------------------------------------------------------*/

static Binan *badbinan(Subarray *sub);

static double *new_calpar(Binan *ban, int oldsize, int newsize);
static double *del_calpar(Binan *ban);
static int fix_calpar(Binan *ban, int nstat, int *t_keep);
static void thr_calpar(Binan *ban, int nstat);

static double *new_orbpar(Binan *ban, int oldsize, int newsize);
static double *del_orbpar(Binan *ban);
static int fix_orbpar(Binan *ban, int nstat, int *t_keep);
static void thr_orbpar(Binan *ban, int nstat);

static Bintel *new_Bintel(Binan *ban, int old_nstat, int new_nstat);
static Bintel *del_Bintel(Binan *ban);
static int fix_Bintel(Binan *ban, int nstat, int *t_keep);

/*.......................................................................
 * Create or re-size a container of binary AIPS AN table info.
 *
 * Input:
 *  sub  Subarray * Pointer to the descriptor of the sub-array that will
 *                  or already does contain the table. Only sub->binan
 *                  and sub->nstat will be used. If NULL, then a new
 *                  container will be allocated here. Otherwise the
 *                  existing container will be re-sized.
 *  nstat     int   The target number of stations to be recorded.
 *  nopcal    int   The required new number of polarization parameters.
 *  numorb    int   The required new number of orbital parameters.
 * Output:
 *  return  Binan * Pointer to the new descriptor, or NULL on error.
 */
Binan *new_Binan(Subarray *sub, int nstat, int nopcal, int numorb)
{
  Binan *ban;    /* The pointer to the new descriptor */
  int old_nstat; /* Previous number of stations handled */
/*
 * Create an initializer for Binan objects.
 */
  static Binan init = {
    0.0, 0.0, 0.0, /* arrayx, arrayy, arrayz */
    0.0,           /* gstia0 */
    0.0,           /* degpdy */
    0.0,           /* freq */
    0.0, 0.0,      /* polarx, polary */
    0.0, 0.0,      /* ut1utc, datutc */
    NULL,          /* calpar */
    NULL,          /* orbpar */
    NULL,          /* bt */
    0, 0,          /* nopcal, numorb */
    "",            /* arrnam[32] */
    "",            /* poltype[16] */
    "",            /* timsys[4] */
    ""             /* rdate[9] */
  };
/*
 * Sanity checks.
 */
  if(sub==NULL) {
    lprintf(stderr, "new_Binan: NULL Subarray descriptor intercepted.\n");
    return NULL;
  };
  if(nstat <= 0) {
    lprintf(stderr, "new_Binan: 0 or -ve number of stations requested.\n");
    return NULL;
  };
/*
 * Get the existing container.
 */
  ban = sub->binan;
/*
 * Enforce sensible values on nopcal and numorb.
 */
  nopcal = nopcal < 0 ? 0 : nopcal;
  numorb = numorb < 0 ? 0 : numorb;
/*
 * How many stations were there previously in the sub-array?
 */
  old_nstat = sub->nstat < 0 ? 0:sub->nstat;
/*
 * New container required?
 */
  if(ban==NULL || old_nstat==0) {
    old_nstat = 0;
    sub->binan = ban = malloc(sizeof(Binan));
    if(ban==NULL)
      return badbinan(sub);
/*
 * Clear members of the new container.
 */
    *ban = init;
  };
/*
 * Allocate/resize and initialize the calpar array.
 */
  if(new_calpar(ban, 2 * ban->nopcal * old_nstat, 2 * nopcal * nstat)==NULL &&
     nopcal != 0)
    return badbinan(sub);
  else
    ban->nopcal = nopcal;
/*
 * Allocate/resize and initialize the orbpar array.
 */
  if(new_orbpar(ban, ban->numorb * old_nstat, numorb * nstat)==NULL &&
     numorb != 0)
    return badbinan(sub);
  else
    ban->numorb = numorb;
/*
 * Allocate/resize and initialize the Bintel *bt array.
 */
  if(new_Bintel(ban, old_nstat, nstat)==NULL)
    return badbinan(sub);
/*
 * Return the revised container.
 */
  return ban;
}

/*.......................................................................
 * Private error reporting and cleanup function of new_Binan().
 */
static Binan *badbinan(Subarray *sub)
{
  lprintf(stderr,"new_Binan: Insufficient memory to record AIPS AN details.\n");
  return del_Binan(sub);
}

/*.......................................................................
 * Delete the Binan descriptor in a given subarray.
 *
 * Input:
 *  sub  Subarray *  Pointer to the descriptor of the sub-array that
 *                   contains the descriptor to be deleted.
 * Output:
 *  return  Binan *  Allways NULL.
 */
Binan *del_Binan(Subarray *sub)
{
  if(sub) {
    Binan *ban = sub->binan;
    if(ban) {
      ban->orbpar = del_orbpar(ban);
      ban->calpar = del_calpar(ban);
      ban->bt = del_Bintel(ban);
      free(sub->binan);
      sub->binan = NULL;
    };
  };
/*
 * Return the deleted descriptor.
 */
  return NULL;
}

/*.......................................................................
 * Remove all but the telescopes flagged as wanted in t_keep[] from
 * a Binan descriptor.
 *
 * Input:
 *  sub     Subarray *  The descriptor of the Subarray containing the
 *                      Binan descriptor to be fixed.
 *  t_keep       int *  An array of nstat elements specifying which
 *                      antennas are to be kept.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int fix_Binan(Subarray *sub, int *t_keep)
{
  Binan *ban;   /* Pointer to the Binan container to be fixed */
  int count=0;  /* The number of telescopes to be kept */
  int i;
/*
 * Sanity checks.
 */
  if(sub==NULL || t_keep==NULL) {
    lprintf(stderr, "fix_Binan: Inapropriate arguments.\n");
    return 1;
  };
/*
 * Get the Binan container.
 */
  ban = sub->binan;
/*
 * Nothing to fix?
 */
  if(ban==NULL)
    return 0;
/*
 * Count the number of stations to be kept.
 */
  for(i=0; i<sub->nstat; i++) {
    if(t_keep[i])
      count++;
  };
/*
 * If there are no stations left, delete the descriptor.
 */
  if(count==0) {
    sub->binan = del_Binan(sub);
    return 0;
  };
/*
 * Fix the orbpar array.
 */
  if(fix_orbpar(ban, sub->nstat, t_keep))
    return 1;
/*
 * Fix the calpar array.
 */
  if(fix_calpar(ban, sub->nstat, t_keep))
    return 1;
/*
 * Fix the Bintel *bt array.
 */
  if(fix_Bintel(ban, sub->nstat, t_keep))
    return 1;
/*
 * Re-thread the orbparm arrays from the ban->orbpar array.
 */
  thr_orbpar(ban, sub->nstat);
/*
 * Re-thread the polcala and polcalb arrays from the ban->calpar array.
 */
  thr_calpar(ban, sub->nstat);
  return 0;
}

/*.......................................................................
 * Private function of new_Binan() used to allocate/resize an array
 * of Bintel AIPS binary AN table telescope descriptors.
 *
 * Input:
 *  ban     Binan * The descriptor that will or does contain a Bintel
 *                  array.
 *  old_nstat int   The number of stations in the existing container.
 *  new_nstat int   The target number of stations to be recorded.
 *  nopcal    int   The required new number of polarization parameters.
 * Output:
 *  return Bintel * The pointer to the revised Bintel array, or NULL on error.
 */
static Bintel *new_Bintel(Binan *ban, int old_nstat, int new_nstat)
{
  Bintel *bt;   /* The array of FITS telescope descriptors */
  Bintel *bp;   /* Pointer into ban->bt[] array */
  int i;
/*
 * Create an initializer for Bintel objects.
 */
  static Bintel init = {
    {0.0,0.0,0.0},                  /* stabxyz[3] */
    NULL,                           /* orbparm[] */
    0.0,                            /* staxof */
    0.0, 0.0,                       /* polaa, polab */
    NULL, NULL,                     /* polcala, polcalb */
    0,                              /* mntsta */
    0,                              /* nosta */
    'R','L',                        /* poltya, poltyb */
    ""                              /* anname[16] */
  };
/*
 * Get a pointer to the array of FITS telescope descriptors.
 */
  bt = ban->bt;
/*
 * Do we require a totally new array of descriptors?
 */
  if(bt==NULL || old_nstat==0) {
    ban->bt = bt = malloc(sizeof(Bintel) * new_nstat);
    if(bt==NULL)
      return bt;
  }
/*
 * Re-size an existing array of descriptors?
 */
  else if(old_nstat != new_nstat) {
    bt = realloc(ban->bt, sizeof(Bintel) * new_nstat);
    if(bt)
      ban->bt = bt;
    else if(new_nstat > old_nstat)
      return NULL;
  };
/*
 * Initialize only the new members of the Bintel array.
 */
  if(new_nstat > old_nstat) {
    bp = &bt[old_nstat];
    for(i=old_nstat; i<new_nstat; i++,bp++)
      *bp = init;
  };
/*
 * The orbparm array of nopcal double's is taken from the
 * ban->orbpar array. Re-thread this array into all the descriptors.
 */
  thr_orbpar(ban, new_nstat);
/*
 * The polcala and polcalb arrays of nopcal double's are taken from the
 * ban->calpar array. Re-thread this array into all the descriptors.
 */
  thr_calpar(ban, new_nstat);
/*
 * Return the revised array.
 */
  return bt;
}

/*.......................................................................
 * Private function of new_Binan() used to delete its contained array
 * of Bintel descriptors. After freeing ban->bt it assigns NULL to it.
 *
 * Input:
 *  ban     Binan *  The container of the Bintel *bt array.
 * Output:
 *  return Bintel *  Allways NULL.
 */
static Bintel *del_Bintel(Binan *ban)
{
  if(ban) {
    if(ban->bt)
      free(ban->bt);
    ban->bt = NULL;
  };
/*
 * Return the deleted descriptor.
 */
  return NULL;
}

/*.......................................................................
 * Private function of fix_Binant() to remove all but the telescopes
 * flagged as wanted in t_keep[].
 *
 * Input:
 *  ban        Binan *  The Binan descriptor containing the Bintel *bt
 *                      array to be fixed.
 *  nstat        int    The existing number of stations.
 *  t_keep       int *  An array of nstat elements specifying which
 *                      antennas are to be kept.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
static int fix_Bintel(Binan *ban, int nstat, int *t_keep)
{
  Bintel *bt;   /* Pointer to the Bintel descriptor in Binan */
/*
 * Get the descriptor to be fixed.
 */
  bt = ban->bt;
/*
 * Anything to be fixed?
 */
  if(bt) {
    int i;
    Bintel *orig;  /* Pointer to original un-fixed element */
    Bintel *dest;  /* Destination for next wanted element */
/*
 * Shuffle the wanted descriptors down the array of decriptors.
 */
    orig = dest = bt;
    for(i=0; i<nstat; i++,orig++) {
      if(t_keep[i])
	*dest++ = *orig;
    };
  };
  return 0;
}

/*.......................................................................
 * Private function of new_Binan() used to allocate/resize a calibration
 * parameter array.
 *
 * Input:
 *  ban     Binan *   The descriptor that will or does contain a calpar
 *                    array.
 *  oldsize   int     The old size for the calpar array.
 *  newsize           The new size for the calpar array.
 * Output:
 *  return double *   The pointer to the calpar array, or NULL on error.
 */
static double *new_calpar(Binan *ban, int oldsize, int newsize)
{
  int i;
/*
 * No existing calpar array?
 */
  if(ban->calpar == NULL)
    oldsize = 0;
/*
 * Calpar array no longer required?
 */
  if(newsize==0) {
    if(ban->calpar)
      free(ban->calpar);
  }
/*
 * New calpar array?
 */
  else if(oldsize==0) {
    ban->calpar = malloc(sizeof(double) * newsize);
    if(ban->calpar == NULL)
      return NULL;
/*
 * Re-size existing calpar array?
 */
  } else if(newsize!=oldsize) {
    double *calpar = realloc(ban->calpar, sizeof(double) * newsize);
    if(calpar)
      ban->calpar = calpar;
    else if(newsize>oldsize)
      return NULL;
  };
/*
 * Initialize just the newly allocated members.
 */
  if(newsize > oldsize) {
    for(i=oldsize; i<newsize; i++)
      ban->calpar[i] = 0.0;
  };
/*
 * Return the revised array.
 */
  return ban->calpar;
}

/*.......................................................................
 * Private function of del_Binan() used to delete its contained calpar
 * array. After freeing ban->calpar it assigns NULL to this member and
 * sets ban->nopcal = 0.
 *
 * Input:
 *  ban     Binan *  The container of the calpar array.
 * Output:
 *  return double *  Allways NULL.
 */
static double *del_calpar(Binan *ban)
{
  if(ban) {
    if(ban->calpar)
      free(ban->calpar);
    ban->calpar = NULL;
    ban->nopcal = 0;
  };
/*
 * Return the deleted calpar array.
 */
  return NULL;
}

/*.......................................................................
 * Private function of fix_Binant() to remove all but the telescopes
 * flagged as wanted in t_keep[].
 *
 * Input:
 *  ban        Binan *  The Binan descriptor containing the calpar
 *                      array to be fixed.
 *  nstat        int    The existing number of stations.
 *  t_keep       int *  An array of nstat elements specifying which
 *                      antennas are to be kept.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
static int fix_calpar(Binan *ban, int nstat, int *t_keep)
{
  double *calpar;   /* Pointer to the array to be fixed */
  int nopcal;       /* The number of cal parameters per feed */
  int i,j;
/*
 * Get a pointer to the calpar array.
 */
  calpar = ban->calpar;
  nopcal = ban->nopcal;
/*
 * Anything to be fixed?
 */
  if(calpar) {
    double *orig;  /* Pointer to original un-fixed entry */
    double *dest;  /* Pointer to destination for next wanted telescope entry */
/*
 * Shuffle the wanted entries down the array. Each station has
 * 2 feeds and nopcal entries per feed.
 */
    orig = dest = calpar;
    for(i=0; i<nstat; i++,orig += 2*nopcal) {
      if(t_keep[i]) {
	for(j=0; j<2*nopcal; j++)
	  *dest++ = *orig;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * The polcala and polcalb arrays of nopcal double's are taken from the
 * ban->calpar array. Re-thread this array into all the ban->bt array
 * of descriptors.
 *
 * Input:
 *  ban    Binan *   The container of the AIPS binary AN table descriptor.
 *  nstat    int     The number of stations recorded in the container.
 */
static void thr_calpar(Binan *ban, int nstat)
{
  if(ban) {
    double *calpar = ban->calpar;
    Bintel *bt = ban->bt;
    int nopcal = ban->nopcal;
    if(bt) {
      int i;
/*
 * Anything to thread?
 */
      if(calpar && nopcal>0) {
	for(i=0; i<nstat; i++,bt++) {
	  bt->polcala = calpar;
	  calpar += nopcal;
	  bt->polcalb = calpar;
	  calpar += nopcal;
	};
/*
 * There is no calpar array to be threaded.
 */
      } else {
	for(i=0; i<nstat; i++,bt++) {
	  bt->polcala = NULL;
	  bt->polcalb = NULL;
	};
      };
    };
  };
  return;
}

/*.......................................................................
 * Private function of new_Binan() used to allocate/resize an orbital
 * parameter array.
 *
 * Input:
 *  ban     Binan *   The descriptor that will or does contain a orbpar
 *                    array.
 *  oldsize   int     The old size for the orbpar array.
 *  newsize           The new size for the orbpar array.
 * Output:
 *  return double *   The pointer to the orbpar array, or NULL on error.
 */
static double *new_orbpar(Binan *ban, int oldsize, int newsize)
{
  int i;
/*
 * No existing orbpar array?
 */
  if(ban->orbpar == NULL)
    oldsize = 0;
/*
 * Orbpar array no longer required?
 */
  if(newsize==0) {
    if(ban->orbpar)
      free(ban->orbpar);
  }
/*
 * New orbpar array?
 */
  else if(oldsize==0) {
    ban->orbpar = malloc(sizeof(double) * newsize);
    if(ban->orbpar == NULL)
      return NULL;
/*
 * Re-size existing orbpar array?
 */
  } else if(newsize!=oldsize) {
    double *orbpar = realloc(ban->orbpar, sizeof(double) * newsize);
    if(orbpar)
      ban->orbpar = orbpar;
    else if(newsize>oldsize)
      return NULL;
  };
/*
 * Initialize just the newly allocated members.
 */
  if(newsize > oldsize) {
    for(i=oldsize; i<newsize; i++)
      ban->orbpar[i] = 0.0;
  };
/*
 * Return the revised array.
 */
  return ban->orbpar;
}

/*.......................................................................
 * Private function of del_Binan() used to delete its contained orbpar
 * array. After freeing ban->orbpar it assigns NULL to this member and
 * sets ban->numorb = 0.
 *
 * Input:
 *  ban     Binan *  The container of the orbpar array.
 * Output:
 *  return double *  Allways NULL.
 */
static double *del_orbpar(Binan *ban)
{
  if(ban) {
    if(ban->orbpar)
      free(ban->orbpar);
    ban->orbpar = NULL;
    ban->numorb = 0;
  };
/*
 * Return the deleted orbpar array.
 */
  return NULL;
}

/*.......................................................................
 * Private function of fix_Binant() to remove the orbital parameter
 * arrays of all but the telescopes that are flagged as wanted in
 * t_keep[].
 *
 * Input:
 *  ban        Binan *  The Binan descriptor containing the orbpar
 *                      array to be fixed.
 *  nstat        int    The existing number of stations.
 *  t_keep       int *  An array of nstat elements specifying which
 *                      antennas are to be kept.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
static int fix_orbpar(Binan *ban, int nstat, int *t_keep)
{
  double *orbpar;   /* Pointer to the array to be fixed */
  int numorb;       /* The number of orbital parameters per station */
  int i,j;
/*
 * Get a pointer to the orbpar array.
 */
  orbpar = ban->orbpar;
  numorb = ban->numorb;
/*
 * Anything to be fixed?
 */
  if(orbpar && numorb>0) {
    double *orig;  /* Pointer to original un-fixed entry */
    double *dest;  /* Pointer to destination for next wanted telescope entry */
/*
 * Shuffle the wanted entries down the array.
 */
    orig = dest = orbpar;
    for(i=0; i<nstat; i++,orig += numorb) {
      if(t_keep[i]) {
	for(j=0; j<numorb; j++)
	  *dest++ = *orig;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * The orbparm arrays of nopcal double's are taken from the
 * ban->orbpar array. Re-thread this array into all the ban->bt array
 * of descriptors.
 *
 * Input:
 *  ban    Binan *   The container of the AIPS binary AN table descriptor.
 *  nstat    int     The number of stations recorded in the container.
 */
static void thr_orbpar(Binan *ban, int nstat)
{
  if(ban) {
    double *orbpar = ban->orbpar;
    Bintel *bt = ban->bt;
    int numorb = ban->numorb;
    if(bt) {
      int i;
/*
 * Anything to thread?
 */
      if(orbpar && numorb>0) {
	for(i=0; i<nstat; i++,bt++) {
	  bt->orbparm = orbpar;
	  orbpar += numorb;
	};
/*
 * There is no orbpar array to be threaded.
 */
      } else {
	for(i=0; i<nstat; i++,bt++)
	  bt->orbparm = NULL;
      };
    };
  };
  return;
}
