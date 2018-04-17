/* FITS primary header HDU code */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "sysfits.h"
#include "fits.h"
#include "utils.h"
#include "phdu.h"
#include "fitkey.h"

static GETFN(get_phdu);
static NEWFN(new_phdu);
static DELFN(del_phdu);
static SIZEFN(siz_phdu);
static ADDFN(add_phdu);
static COPYFN(cop_phdu);
static ENDFN(end_phdu);

Hdutab phdufns={"IMAGE", "IMAGE", get_phdu, new_phdu, del_phdu, siz_phdu,
		add_phdu, cop_phdu, end_phdu};

static long imagesize(Phdu *phdu);

static Hdu *newimage(Hdutype type, Bitpix bitpix, int naxis, int *dims,
		     int groups, int pcount, int gcount, char *extname,
		     int extver, int extlevel);

static Gpar *new_Gpars(Phdu *phdu);
static Imaxis *new_axes(Phdu *phdu);

/*
 * Declare a scratch buffer used to collect and send the different scale and
 * offset factors of different random group parameters to arrconv().
 * It is independantly used by wgroup() and rgroup().
 */
enum {PHDU_NPAR=30};
static Offscal offscal[PHDU_NPAR];

/*.......................................................................
 * Clear the derived part of a primary header descriptor to the
 * point at which it can be safely sent to del_Hdu().
 *
 * Input:
 *  hdu     Hdu *  The base-class pointer to a Phdu descriptor.
 * Output:
 *  return  Hdu *  'hdu', or NULL if an error occurs and hdu is deleted.
 */
static NEWFN(new_phdu)
{
  Phdu *phdu = (Phdu *) hdu;
  phdu->origin = NULL;
  phdu->date_obs = NULL;
  phdu->telescop = NULL;
  phdu->instrume = NULL;
  phdu->observer = NULL;
  phdu->object = NULL;
  phdu->author = NULL;
  phdu->referenc = NULL;
  phdu->bunit  = NULL;
  phdu->axes   = NULL;
  phdu->pars   = NULL;
  phdu->bscale = 1.0;
  phdu->bzero = 0.0;
  phdu->blank = NONULL;
  phdu->equinox = 0.0;
  phdu->datamin = phdu->datamax = 0.0;
  return hdu;
}

/*.......................................................................
 * Delete the derived parts of a primary or image HDU descriptor.
 *
 * Input:
 *  hdu      Hdu *  The base-class pointer to the Phdu descriptor.
 */
static DELFN(del_phdu)
{
  Phdu *phdu = (Phdu *) hdu;
  Imaxis *axis;
  Gpar *gpar;
  int i;
  if(phdu->origin)
    free(phdu->origin);
  if(phdu->date_obs)
    free(phdu->date_obs);
  if(phdu->telescop)
    free(phdu->telescop);
  if(phdu->instrume)
    free(phdu->instrume);
  if(phdu->observer)
    free(phdu->observer);
  if(phdu->object)
    free(phdu->object);
  if(phdu->author)
    free(phdu->author);
  if(phdu->referenc)
    free(phdu->referenc);
  if(phdu->bunit)
    free(phdu->bunit);
/*
 * Clean-up the array of axis descriptors.
 */
  if(phdu->axes != 0) {
    for(i=0; i<phdu->naxis; i++) {
      axis = &phdu->axes[i];
      if(axis->ctype)
	free(axis->ctype);
    };
    free(phdu->axes);
  };
/*
 * CLean-up the array of group-parameter descriptors.
 */
  if(phdu->pars != 0) {
    for(i=0; i<phdu->pcount; i++) {
      gpar = &phdu->pars[i];
      if(gpar->ptype)
	free(gpar->ptype);
    };
    free(phdu->pars);
  };
  return;
}

/*.......................................................................
 * Read the PRIMARY header or IMAGE extension specific header keywords
 * and record them in a given primary-HDU descriptor.
 *
 * Input:
 *  fits  Fits *  The descriptor of the FITS file.
 * Input/Output:
 *  hdu    Hdu *  Base-class pointer to a Phdu descriptor.
 * Output:
 *  return int    0 - OK.
 */
static GETFN(get_phdu)
{
  Phdu *phdu = (Phdu *) hdu;
  Fitkey key;        /* Keyword-value pair temporary storage */
  int first=1;       /* False after first time through a loop */
  int i;
  /* Axis descriptor keywords */
  enum {PHDU_CTYPE, PHDU_CRPIX, PHDU_CRVAL, PHDU_CDELT, PHDU_CROTA};
  Fitkey axkeys[]={
    {"CTYPE", 0, PHDU_CTYPE, DAT_STR, NULL, NULL},
    {"CRPIX", 0, PHDU_CRPIX, DAT_DBL, NULL, NULL},
    {"CRVAL", 0, PHDU_CRVAL, DAT_DBL, NULL, NULL},
    {"CDELT", 0, PHDU_CDELT, DAT_DBL, NULL, NULL},
    {"CROTA", 0, PHDU_CROTA, DAT_DBL, NULL, NULL}
  };
  /* Random-groups parameter descriptor keywords */
  enum {PHDU_PTYPE, PHDU_PSCAL, PHDU_PZERO};
  Fitkey pkeys[]={
    {"PTYPE", 0, PHDU_PTYPE, DAT_STR, NULL, NULL},
    {"PSCAL", 0, PHDU_PSCAL, DAT_DBL, NULL, NULL},
    {"PZERO", 0, PHDU_PZERO, DAT_DBL, NULL, NULL}
  };
  /* Optional keywords */
  enum {PHDU_BSCALE, PHDU_BZERO, PHDU_BUNIT, PHDU_BLANK, PHDU_OBJECT,
	PHDU_TELESCOP, PHDU_ORIGIN, PHDU_DATE_OBS, PHDU_INSTRUME,
	PHDU_OBSERVER, PHDU_AUTHOR, PHDU_REFERENC, PHDU_EQUINOX, PHDU_EPOCH,
	PHDU_DATAMAX, PHDU_DATAMIN};
  Fitkey okeys[]={
    {"BSCALE",   0, PHDU_BSCALE,   DAT_DBL, NULL, NULL},
    {"BZERO",    0, PHDU_BZERO,    DAT_DBL, NULL, NULL},
    {"BUNIT",    0, PHDU_BUNIT,    DAT_STR, NULL, NULL},
    {"BLANK",    0, PHDU_BLANK,    DAT_INT, NULL, NULL},
    {"OBJECT",   0, PHDU_OBJECT,   DAT_STR, NULL, NULL},
    {"TELESCOP", 0, PHDU_TELESCOP, DAT_STR, NULL, NULL},
    {"ORIGIN",   0, PHDU_ORIGIN,   DAT_STR, NULL, NULL},
    {"DATE-OBS", 0, PHDU_DATE_OBS, DAT_STR, NULL, NULL},
    {"INSTRUME", 0, PHDU_INSTRUME, DAT_STR, NULL, NULL},
    {"OBSERVER", 0, PHDU_OBSERVER, DAT_STR, NULL, NULL},
    {"AUTHOR",   0, PHDU_AUTHOR,   DAT_STR, NULL, NULL},
    {"REFERENC", 0, PHDU_REFERENC, DAT_STR, NULL, NULL},
    {"EQUINOX",  0, PHDU_EQUINOX,  DAT_DBL, NULL, NULL},
    {"EPOCH",    0, PHDU_EPOCH,    DAT_DBL, NULL, NULL},
    {"DATAMAX",  0, PHDU_DATAMAX,  DAT_DBL, NULL, NULL},
    {"DATAMIN",  0, PHDU_DATAMIN,  DAT_DBL, NULL, NULL},
  };
/*
 * Check the input Hdu type.
 */
  if(hdu->type != F_IMAGE && hdu->type != F_PRIMARY) {
    fprintf(stderr, "get_phdu: syserror: Incompatible HDU descriptor received\n");
    return 1;
  };
/*
 * Get the optional keywords.
 */
  while(next_key(fits, hdu, okeys, sizeof(okeys)/sizeof(Fitkey), EOH_SEEK,
		 &key) == 0) {
    switch(key.keyid) {
    case PHDU_BSCALE:
      phdu->bscale = KEYDBL(key);
      break;
    case PHDU_BZERO:
      phdu->bzero = KEYDBL(key);
      break;
    case PHDU_BUNIT:
      phdu->bunit = fitsstr(KEYSTR(key));
      break;
    case PHDU_BLANK:
      phdu->blank = KEYINT(key);
      break;
    case PHDU_OBJECT:
      phdu->object = fitsstr(KEYSTR(key));
      break;
    case PHDU_TELESCOP:
      phdu->telescop = fitsstr(KEYSTR(key));
      break;
    case PHDU_ORIGIN:
      phdu->origin = fitsstr(KEYSTR(key));
      break;
    case PHDU_DATE_OBS:
      phdu->date_obs = fitsstr(KEYSTR(key));
      break;
    case PHDU_INSTRUME:
      phdu->instrume = fitsstr(KEYSTR(key));
      break;
    case PHDU_OBSERVER:
      phdu->observer = fitsstr(KEYSTR(key));
      break;
    case PHDU_AUTHOR:
      phdu->author = fitsstr(KEYSTR(key));
      break;
    case PHDU_REFERENC:
      phdu->referenc = fitsstr(KEYSTR(key));
      break;
    case PHDU_EQUINOX: case PHDU_EPOCH:
      phdu->equinox = KEYDBL(key);
      break;
    case PHDU_DATAMAX:
      phdu->datamax = KEYDBL(key);
      break;
    case PHDU_DATAMIN:
      phdu->datamin = KEYDBL(key);
      break;
    default:
      break;
    };
  };
/*
 * Allocate an array of NAXIS axis descriptors.
 */
  if(phdu->naxis > 0 && !phdu->axes && !new_axes(phdu))
    return 1;
/*
 * See if there are any axis descriptor keywords. If so allocate an
 * array of NAXIS axis descriptors and fill it with the contents
 * of the header.
 */
  if(phdu->naxis > 0) {
    first=1;
    hdu->nextline = 0;
    while(next_key(fits, hdu, axkeys, sizeof(axkeys)/sizeof(Fitkey),
		   EOH_SEEK, &key) == 0) {
/*
 * Install the keyword value in the appropriate slot.
 */
      if(key.extn > 0 && key.extn <= phdu->naxis) {
	i = key.extn - 1;
	switch(key.keyid) {
	case PHDU_CTYPE:
	  phdu->axes[i].ctype = fitsstr(KEYSTR(key));
	  break;
	case PHDU_CRPIX:
	  phdu->axes[i].crpix = KEYDBL(key);
	  break;
	case PHDU_CRVAL:
	  phdu->axes[i].crval = KEYDBL(key);
	  break;
	case PHDU_CDELT:
	  phdu->axes[i].cdelt = KEYDBL(key);
	  break;
	case PHDU_CROTA:
	  phdu->axes[i].crota = KEYDBL(key);
	  break;
	};
      };
    };
  };
/*
 * See if there are any group parameter descriptor keywords. If so allocate
 * an array of pcount axis descriptors and fill it with the contents
 * of the header.
 */
  if(phdu->groups && phdu->pcount > 0) {
    first=1;
    hdu->nextline = 0;
    while(next_key(fits, hdu, pkeys, sizeof(pkeys)/sizeof(Fitkey), EOH_SEEK,
		   &key) == 0)
    {
/*
 * Allocate the array of descriptors if necessary.
 */
      if(!phdu->pars && new_Gpars(phdu)==NULL)
	break;
/*
 * Install the keyword value in the appropriate slot.
 */
      if(key.extn > 0 && key.extn <= phdu->pcount) {
	i = key.extn - 1;
	switch(key.keyid) {
	case PHDU_PTYPE:
	  phdu->pars[i].ptype = fitsstr(KEYSTR(key));
	  break;
	case PHDU_PSCAL:
	  phdu->pars[i].pscal = KEYDBL(key);
	  break;
	case PHDU_PZERO:
	  phdu->pars[i].pzero = KEYDBL(key);
	  break;
	};
      };
    };
  };
/*
 * Record the number of objects in the image-array of each group.
 */
  phdu->imsize = imagesize(phdu);
  return 0;
}

/*.......................................................................
 * Return the size of an IMAGE HDU descriptor.
 *
 * Output:
 *  return  size_t  sizeof(Phdu).
 */
static SIZEFN(siz_phdu)
{
  return sizeof(Phdu);
}

/*.......................................................................
 * Find the HDU of an IMAGE extension by its name and version number.
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
 *  return   Phdu *  The pointer to the IMAGE HDU descriptor
 *                   or NULL if not found.
 */
Phdu *find_image(Fits *fits, char *extname, int extver, Hdu *prev)
{
  return (Phdu *) find_hdu(fits, F_IMAGE, extname, extver, prev);
}

/*.......................................................................
 * Read a given random groups entry from an IMAGE HDU.
 *
 * Input:
 *  fits      Fits *  The FITS file descriptor.
 *  phdu      Phdu *  The descriptor of the IMAGE HDU to be read from.
 *  igroup    long    The 0-relative index of the group to be read from.
 *  start     long    The 0-relative index of the first parameter to be
 *                    read from the group.
 *  nobj      long    The maximum number of parameters to read.
 *  type   Fittype    The declared type of array data[] - conversion to this
 *                    type will be performed if meaningful.
 *  doscale    int    If true apply offset and scale factors.
 *  flags Fitsflag *  If not NULL, then this should point to an array
 *                    of nobj elements to be used to record which elements
 *                    of data[] are blanked 1=blanked, 0=OK.
 * Output:
 *  data      void *  Send an array of type 'type' and with a dimension of
 *                    at least 'nobj'. This will contain the output data.
 *  return    long    The number of objects read.
 */
long rgroup(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	       Fittype type, int doscale, Fitsflag *flags, void *data)
{
  Hdu *hdu = (Hdu *) phdu;
  long offset;/* Offset into data-segment of HDU (FITS bytes) */
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL || fits==NULL || data==NULL) {
    fprintf(stderr, "rgroup: NULL parameter intercepted\n");
    return 0L;
  };
/*
 * Is this a random-groups HDU?
 */
  if(!hdu->groups) {
    fprintf(stderr, "rgroup: The given HDU does not contain random-groups.\n");
    return 0L;
  };
/*
 * The group index may exceed the currently recorded number of groups
 * in the file only if this is a random-groups HDU whose data segment
 * is incomplete.
 */
  if(igroup < 0 || (hdu->state!=HDU_DATA && igroup >= hdu->gcount)) {
    fprintf(stderr, "rgroup: Group index (%ld) out of range.\n", igroup+1L);
    return 0L;
  };
  if(start < 0 || start >= hdu->pcount) {
    fprintf(stderr, "rgroup: Group-paremeter start index out of range\n");
    return 0L;
  };
/*
 * Does this HDU contain any random groups?
 */
  if(hdu->pcount==0)
    return 0L;
/*
 * Compute the offset to the start of the group-parameter block wrt the
 * start of the data segment.
 */
  offset = hdu->grpsize * igroup;
/*
 * Prevent field overflow.
 */
  if(nobj+start > hdu->pcount)
    nobj = hdu->pcount-start;
/*
 * Different random-group parameters have different offset and scale
 * factors. If scaling has been requested, assemble an array of these
 * to send to get_data().
 */
  if(doscale) {
    Gpar *par;   /* Descriptor of group parameter to be processed next */
    size_t size; /* The number of bytes in an element of data[] */
    int obj;     /* The index of the object being read wrt 'start' */
    int i;
/*
 * Get the number of bytes per element of data[].
 */
    size = machsize(type);
/*
 * Read the data in chunks of up to PHDU_NPAR elements.
 */
    par = &phdu->pars[start];
    for(obj=0; obj<nobj; obj+=PHDU_NPAR) {
/*
 * How many elements can be processed in this pass?
 */
      int nnew = (nobj-obj < PHDU_NPAR) ? nobj-obj : PHDU_NPAR;
/*
 * Assemble a list of the next 'nnew' random-group parameter scale and
 * offset factors.
 */
      for(i=0; i<nnew; i++,par++) {
	offscal[i].off = par->pzero;
	offscal[i].mul = par->pscal;
      };
/*
 * Get the next 'nnew' parameters.
 */
      if(get_data(fits, hdu, offset, dat_type(hdu), start+obj, nnew, type,
		  0.0, 1.0, offscal, flags, phdu->blank,
		  ((char *)data) + obj*size))
	return 0L;
    };
  }
/*
 * No array of scale factors is required if scaling was not requested.
 */
  else {
    if(get_data(fits, hdu, offset, dat_type(hdu), start, nobj, type,
		0.0, 1.0, NULL, flags, phdu->blank, data))
      return 0L;
  };
  return nobj;
}

/*.......................................................................
 * Write a given random groups entry to an IMAGE HDU.
 *
 * Input:
 *  fits      Fits * The FITS file descriptor.
 *  phdu      Phdu * The descriptor of the IMAGE HDU to be written into.
 *  igroup    long   The 0-relative index of the group to be written into.
 *  start     long   The 0-relative index of the first parameter to be
 *                   written from the group.
 *  nobj      long   The maximum number of parameters to write.
 *  type   Fittype   The declared type of array data[] - conversion from this
 *                   type will be performed if meaningful.
 *  doscale    int   If true remove offset and scale factors.
 *  flags Fitsflag * If not NULL, then this should point to an array
 *                   of nobj elements specifying which elements
 *                   of data[] are to be blanked 1=blanked, 0=OK.
 *  data      void * Send an array of type 'type' and with a dimension of
 *                   at least 'nobj'. This will contain the output data.
 * Output:
 *  return    long   The number of objects written.
 */
long wgroup(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	    Fittype type, int doscale, Fitsflag *flags, void *data)
{
  Hdu *hdu = (Hdu *) phdu;
  long offset;  /* Offset into data-segment of HDU (FITS bytes) */
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL || fits==NULL || data==NULL) {
    fprintf(stderr, "wgroup: NULL parameter intercepted\n");
    return 0L;
  };
/*
 * Is this a random-groups HDU?
 */
  if(!hdu->groups) {
    fprintf(stderr, "wgroup: The given HDU does not contain random-groups.\n");
    return 0L;
  };
  if(igroup < 0) {
    fprintf(stderr, "wgroup: Negative group index intercepted.\n");
    return 0L;
  };
  if(start < 0 || start >= hdu->pcount) {
    fprintf(stderr, "wgroup: Group-paremeter start index out of range\n");
    return 0L;
  };
/*
 * Is there room for any random group-parameters?
 */
  if(hdu->pcount==0)
    return 0L;
/*
 * Expand the established range of groups?
 */
  if(igroup >= hdu->gcount) {
    if(hdu->state==HDU_DATA) {
      hdu->gcount = igroup+1;
    } else {
      fprintf(stderr, "wgroup: Can't expand GCOUNT to incorporate group %ld.\n",
	      igroup+1L);
      return 0L;
    };
  };
/*
 * Compute the offset to the start of the group-parameter block wrt the
 * start of the data segment.
 */
  offset = hdu->grpsize * igroup;
/*
 * Prevent field overflow.
 */
  if(nobj+start > hdu->pcount)
    nobj = hdu->pcount-start;
/*
 * Different random-group parameters have different offset and scale
 * factors. If scaling has been requested, assemble an array of these
 * to send to put_data().
 */
  if(doscale) {
    Gpar *par;   /* Descriptor of group parameter to be processed next */
    size_t size; /* The number of bytes in an element of data[] */
    int obj;     /* The index of the object being written wrt 'start' */
    int i;
/*
 * Get the number of bytes per element of data[].
 */
    size = machsize(type);
/*
 * Write the data in chunks of up to PHDU_NPAR elements.
 */
    par = &phdu->pars[start];
    for(obj=0; obj<nobj; obj+=PHDU_NPAR) {
/*
 * How many elements can be processed in this pass?
 */
      int nnew = (nobj-obj < PHDU_NPAR) ? nobj-obj : PHDU_NPAR;
/*
 * Assemble a list of the next 'nnew' random-group parameter scale and
 * offset factors.
 */
      for(i=0; i<nnew; i++,par++) {
	offscal[i].off = -par->pzero;
	offscal[i].mul = 1.0/par->pscal;
      };
/*
 * Write the next 'nnew' parameters.
 */
      if(put_data(fits, hdu, offset, dat_type(hdu), start+obj, nnew, type,
		  0.0, 1.0, offscal, flags, phdu->blank,
		  ((char *)data) + obj*size))
	return 0;
    };
  }
/*
 * No array of scale factors is required if scaling was not requested.
 */
  else {
    if(put_data(fits, hdu, offset, dat_type(hdu), start, nobj, type,
		0.0, 1.0, NULL, flags, phdu->blank, data))
      return 0;
  };
  return nobj;
}

/*.......................................................................
 * Read the image-array of a given group of an IMAGE HDU.
 *
 * Input:
 *  fits      Fits *  The FITS file descriptor.
 *  phdu      Phdu *  The descriptor of the IMAGE HDU to be read.
 *  igroup    long    The 0-relative index of the group to be read.
 *                    If this is not a random-groups header send 0L.
 *  start     long    The 0-relative index of the first element to be
 *                    read from the image array.
 *  nobj      long    The maximum number of elements to read.
 *  type   Fittype    The delcared type of array data[] - conversion to this
 *                    type will be performed if meaningful.
 *  doscale    int    If true apply offset and scale factors.
 *  flags Fitsflag *  If not NULL, then this should point to an array
 *                    of nobj elements to be used to record which elements
 *                    of data[] are blanked 1=blanked, 0=OK.
 * Output:
 *  data    void *  Send an array of type 'type' and with a dimension of
 *                  at least 'nobj'. This will contain the output data.
 *  return  long    The number of objects read.
 */
long rimage(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	       Fittype type, int doscale, Fitsflag *flags, void *data)
{
  Hdu *hdu = (Hdu *) phdu;
  long offset;/* Offset into data-segment of HDU (FITS bytes) */
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL || fits==NULL || data==NULL) {
    fprintf(stderr, "rimage: NULL parameter intercepted\n");
    return 0L;
  };
/*
 * The group index may exceed the currently recorded number of groups
 * in the file only if this is a random-groups HDU whose data segment
 * is incomplete.
 */
  if(igroup < 0 ||
     ((hdu->state!=HDU_DATA || !hdu->groups) && igroup >= hdu->gcount)) {
    fprintf(stderr, "rimage: Group index (%ld) out of range.\n", igroup+1L);
    return 0L;
  };
  if(start < 0 || start >= phdu->imsize) {
    fprintf(stderr, "rimage: Image-array start index out of range\n");
    return 0L;
  };
/*
 * Compute the offset to the start of the group wrt the start of the data
 * segment. 
 */
  offset = hdu->grpsize * igroup;
/*
 * Enforce limits on the number of elements that can be read.
 */
  if(nobj>phdu->imsize)
    nobj = phdu->imsize;
/*
 * Get the data.
 */
  if(get_data(fits, hdu, offset, dat_type(hdu), hdu->pcount+start, nobj, type,
	      doscale ? phdu->bzero:0.0, doscale ? phdu->bscale:1.0, NULL,
	      flags, phdu->blank, data))
    return 0;
  return nobj;
}

/*.......................................................................
 * Write the image-array of a given group of an IMAGE HDU.
 *
 * Input:
 *  fits      Fits * The FITS file descriptor.
 *  phdu      Phdu * The descriptor of the IMAGE HDU to be written.
 *  igroup    long   The 0-relative index of the group to be written.
 *                   If this is not a random-groups header send 0L.
 *  start     long   The 0-relative index of the first element to be
 *                   written from the image array.
 *  nobj      long   The maximum number of elements to write.
 *  type   Fittype   The declared type of array data[] - conversion from
 *                   this type will be performed if meaningful.
 *  doscale    int   If true remove offset and scale factors.
 *  flags Fitsflag * If not NULL, then this should point to an array
 *                   of nobj elements specifying which elements
 *                   of data[] are to be blanked 1=blanked, 0=OK.
 *  data      void * Send an array of type 'type' and with a dimension of
 *                   at least 'nobj', containing the data to be written.
 * Output:
 *  return    long   The number of objects written.
 */
long wimage(Fits *fits, Phdu *phdu, long igroup, long start, long nobj,
	       Fittype type, int doscale, Fitsflag *flags, void *data)
{
  Hdu *hdu = (Hdu *) phdu;
  long offset; /* Offset into data-segment of HDU (FITS bytes) */
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL || fits==NULL || data==NULL) {
    fprintf(stderr, "wimage: NULL parameter intercepted\n");
    return 0L;
  };
  if(igroup < 0) {
    fprintf(stderr, "wimage: Negative group index intercepted.\n");
    return 0L;
  };
  if(start < 0 || start >= phdu->imsize) {
    fprintf(stderr, "wimage: Image-array start index out of range\n");
    return 0L;
  };
/*
 * Expand the established range of groups?
 */
  if(igroup >= hdu->gcount) {
    if(hdu->state==HDU_DATA) {
      hdu->gcount = igroup+1;
    } else {
      fprintf(stderr, "wimage: Can't expand GCOUNT to incorporate group %ld.\n",
	      igroup+1L);
      return 0L;
    };
  };
/*
 * Compute the offset to the start of the group wrt the start of the data
 * segment. 
 */
  offset = hdu->grpsize * igroup;
/*
 * Enforce limits on the number of elements that can be written.
 */
  if(nobj+start>phdu->imsize)
    nobj = phdu->imsize-start;
/*
 * Write the data.
 */
  if(put_data(fits, hdu, offset, dat_type(hdu), hdu->pcount+start, nobj, type,
	      doscale ? -phdu->bzero:0.0, doscale ? 1.0/phdu->bscale:1.0,
	      NULL, flags, phdu->blank, data))
    return 0;
  return nobj;
}

/*.......................................................................
 * Return the number of elements in an image-array.
 *
 * Input:
 *  phdu    Phdu *   The descriptor of the image HDU.
 * Output:
 *  return  long     The number of elements in the array.
 */
static long imagesize(Phdu *phdu)
{
  long nobj=1L;
  int i;
/*
 * Determine the number of objects in one group image-array.
 */
  if(phdu->naxis==0 || phdu->dims == NULL) {
    nobj = 0L;
  } else if(phdu->groups) {
    for(i=1; i<phdu->naxis; i++)
      nobj *= phdu->dims[i];
  } else {
    for(i=0; i<phdu->naxis; i++)
      nobj *= phdu->dims[i];
  };
/*
 * Return the number of elements.
 */
  return nobj;
}

/*.......................................................................
 * Private function for use by new_image and new_primary. Create and
 * return a new image hdu.
 *
 * Input:
 *  type    Hdutype   F_PRIMARY or F_IMAGE only.
 *  bitpix   Bitpix   The data-type of the image - see hdu.h for the list
 *                    of valid types.
 *  naxis       int   The number of dimensions to the image array.
 *  dims        int * An array of size 'int dims[naxis]' containing the
 *                    number of elements along each dimension of the
 *                    image array.
 *  groups      int   True if the image is to have random-groups.
 *  pcount      int   The number of misc parameters per group.
 *  gcount      int   The number of groups.
 *  extname    char * If this is to be an IMAGE extension, send a name to
 *                    assign to it. If this is a primary HDU send "PRIMARY".
 *                    A copy of the string will be made.
 *  extver      int   The version number to assign to the extension. If
 *                    0 - then when the HDU is added to a FITS file, it
 *                    will be assigned the next highest version number for
 *                    name 'extname' and level 'extlevel'.
 *  extlevel    int   The hierachical level of this extension. Send 1 or
 *                    0 if not relevant.
 * Output:
 *  return      Hdu * The new HDU descriptor.
 */
static Hdu *newimage(Hdutype type, Bitpix bitpix, int naxis, int *dims,
		     int groups, int pcount, int gcount, char *extname,
		     int extver, int extlevel)
{
  Phdu *phdu;
  Hdu *hdu;
/*
 * Sanity check arguments.
 */
  if(type!=F_PRIMARY && type!=F_IMAGE) {
    fprintf(stderr, "newimage: Illegal HDU type intercepted\n");
    return NULL;
  };
  if(dims==NULL) {
    fprintf(stderr, "newimage: NULL dims argument intercepted\n");
    return NULL;
  };
  if(groups && dims[0]!=0) {
    fprintf(stderr, "newimage: NAXIS1 must be 0 in a random-groups HDU\n");
    return NULL;
  };
/*
 * Allocate and default initialize a Phdu descriptor.
 */
  hdu = new_Hdu(type);
  if(hdu==NULL)
    return hdu;
/*
 * Initialize the base-class part of the HDU descriptor.
 */
  hdu = ini_Hdu(hdu, bitpix, dims, naxis, groups, groups?pcount:0,
		groups?gcount:1, extname, extver, extlevel, 0, 0);
  if(hdu==NULL)
    return hdu;
/*
 * Get the derived version of the descriptor.
 */
  phdu = (Phdu *) hdu;
/*
 * Record the number of objects in the image-array of each group.
 */
  phdu->imsize = imagesize(phdu);
/*
 * Allocate an array of NAXIS axis descriptors.
 */
  if(naxis > 0 && !new_axes(phdu))
    return del_Hdu(hdu);
  return hdu;
}

/*.......................................................................
 * Create and return a new image-extension hdu.
 *
 * Input:
 *  bitpix   Bitpix   The data-type of the image - see hdu.h for the list
 *                    of valid types.
 *  naxis       int   The number of dimensions to the image array.
 *  dims        int * An array of size 'int dims[naxis]' containing the
 *                    number of elements along each dimension of the
 *                    image array.
 *  extname    char * If this is to be an IMAGE extension, send a name to
 *                    assign to it. A copy of the string will be made.
 *  extver      int   The version number to assign to the extension. If
 *                    0 - then when the HDU is added to a FITS file, it
 *                    will be assigned the next highest version number for
 *                    name 'extname' and level 'extlevel'.
 *  extlevel    int   The hierachical level of this extension. Send 1 or
 *                    0 if not relevant.
 * Output:
 *  return      Hdu * The new HDU descriptor.
 */
Hdu *new_image(Bitpix bitpix, int naxis, int *dims, char *extname, int extver,
	       int extlevel)
{
  return newimage(F_IMAGE, bitpix, naxis, dims, 0, 0, 1, extname, extver,
		  extlevel);
}

/*.......................................................................
 * Create and return a new primary HDU.
 *
 * Input:
 *  bitpix   Bitpix   The data-type of the image - see hdu.h for the list
 *                    of valid types.
 *  naxis       int   The number of dimensions to the image array.
 *  dims        int * An array of size 'int dims[naxis]' containing the
 *                    number of elements along each dimension of the
 *                    image array.
 *  groups      int   True if the image is to have random-groups.
 *  pcount      int   The number of misc parameters per group.
 *  gcount      int   The number of groups.
 * Output:
 *  return      Hdu * The new HDU descriptor.
 */
Hdu *new_primary(Bitpix bitpix, int naxis, int *dims, int groups, int pcount,
		 int gcount)
{
  return newimage(F_PRIMARY, bitpix, naxis, dims, groups, pcount, gcount,
		  "PRIMARY", 0, 1);
}

/*.......................................................................
 * Set the description keywords of an axis in a primary or image HDU.
 *
 * Input:
 *  hdu      Hdu *  The HDU descriptor of the primary or image HDU.
 *  axis     int    The 1-relative index of the axis whose description
 *                  is to be assigned.
 *  ctype   char *  A name to refer to the axis by.
 *  crpix double   The reference pixel.
 *  crval double   The value on this axis at pixel 'crpix'.
 *  cdelt double   The value increment per pixel on this axis.
 *  crota double   The rotation angle of the axis - see the NOST standard
 *                 for details.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int setaxis(Hdu *hdu, int axis, char *ctype, double crpix, double crval,
	    double cdelt, double crota)
{
  Phdu *phdu = (Phdu *) hdu;
  Imaxis *faxis;
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL) {
    fprintf(stderr, "setaxis: Intercepted NULL argument HDU descriptor\n");
    return 1;
  };
  if(hdu->complete) {
    fprintf(stderr,
	    "setaxis: Illegal attempt to change established HDU structure\n");
    return 1;
  };
  if(hdu->type != F_IMAGE && hdu->type != F_PRIMARY) {
    fprintf(stderr, "setaxis: Inappropriate HDU type intercepted\n");
    return 1;
  };
  if(axis < 1 || axis > hdu->naxis) {
    fprintf(stderr, "setaxis: Axis index out of range 1->%d\n", hdu->naxis);
    return 1;
  };
/*
 * Assign the axis description values to their slots.
 */
  faxis = &phdu->axes[axis-1];
  if(!faxis->ctype && ctype)
    faxis->ctype = fitsstr(ctype);
  faxis->crpix = crpix;
  faxis->crval = crval;
  faxis->cdelt = cdelt;
  faxis->crota = crota;
  return 0;
}

/*.......................................................................
 * Set the description keywords of a random-groups parameter in a primary
 * or image HDU.
 *
 * Input:
 *  hdu      Hdu *  The HDU descriptor of the primary or image HDU.
 *  ipar     int    The 1-relative index of the parameter whose description
 *                  is to be assigned.
 *  ptype   char *  A name to refer to the parameter by.
 *  pscal double    The scale factor for the parameter.
 *  pzero double    The zero offset of the parameter.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int setgroup(Hdu *hdu, int ipar, char *ptype, double pscal, double pzero)
{
  Phdu *phdu = (Phdu *) hdu;
  Gpar *gpar;
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL) {
    fprintf(stderr, "setgroup: Intercepted NULL argument HDU descriptor\n");
    return 1;
  };
  if(hdu->complete) {
    fprintf(stderr,
	    "setgroup: Illegal attempt to change established HDU structure\n");
    return 1;
  };
  if(hdu->type != F_IMAGE && hdu->type != F_PRIMARY) {
    fprintf(stderr, "setgroup: Inappropriate HDU type intercepted\n");
    return 1;
  };
  if(ipar < 1 || ipar > hdu->pcount) {
    fprintf(stderr, "setgroup: Parameter index out of range 1->%d\n",
	    hdu->pcount);
    return 1;
  };
/*
 * Allocate the array of parameter descriptors if necessary.
 */
  if(!phdu->pars && new_Gpars(phdu)==NULL)
    return 1;
/*
 * Assign the parameter description values to their slots.
 */
  gpar = &phdu->pars[ipar-1];
  if(!gpar->ptype && ptype)
    gpar->ptype = fitsstr(ptype);
  gpar->pscal = pscal;
  gpar->pzero = pzero;
  return 0;
}

/*.......................................................................
 * Set the description keywords of the data array in a primary or image HDU.
 *
 * Input:
 *  hdu        Hdu *  The HDU descriptor of the primary or image HDU.
 *  bscale  double    The scale factor of the data.
 *  bzero   double    The zero-offset of the data.
 *  bunit     char *  The name of the units describing the data.
 *  blank     long    The integer used to represent blanked data values.
 *                    Send blank=NONULL if blanking not required.
 *  datamin double    The minimum data value in the data array.
 *  datamax double    The maximum data value in the data array.
 * Output:
 *  return     int    0 - OK.
 *                   1 - Error.
 */
int setimage(Hdu *hdu, double bscale, double bzero, char *bunit, long blank,
	     double datamin, double datamax)
{
  Phdu *phdu = (Phdu *) hdu;
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL) {
    fprintf(stderr, "setimage: Intercepted NULL argument HDU descriptor\n");
    return 1;
  };
  if(hdu->complete) {
    fprintf(stderr,
	    "setimage: Illegal attempt to change established HDU structure\n");
    return 1;
  };
  if(hdu->type != F_IMAGE && hdu->type != F_PRIMARY) {
    fprintf(stderr, "setimage: Inappropriate HDU type intercepted\n");
    return 1;
  };
/*
 * Assign the description values to their slots.
 */
  phdu->bscale = bscale;
  phdu->bzero = bzero;
  if(!phdu->bunit && bunit)
    phdu->bunit = fitsstr(bunit);
  phdu->blank = blank;
  phdu->datamin = datamin<datamax ? datamin : datamax;
  phdu->datamax = datamin<datamax ? datamax : datamin;
  return 0;
}

/*.......................................................................
 * Allocate an array of random-groups descriptors in a primary or image
 * HDU.
 *
 * Input:
 *  phdu   Phdu *  The PRIMARY or IMAGE extension phdu to allocate the
 *                 array for. phdu->pcount must have been set before
 *                 the call.
 * Output:
 *  return Gpar *  The pointer to the allocated array of descriptors
 *                 or NULL on error. This will also be assigned to
 *                 phdu->pars.
 */
static Gpar *new_Gpars(Phdu *phdu)
{
  int i;
/*
 * Sanity checks.
 */
  if(phdu==NULL) {
    fprintf(stderr, "new_Gpars: NULL Phdu descriptor intercepted\n");
    return NULL;
  };
  if(phdu->type != F_PRIMARY) {
    fprintf(stderr, "new_Gpars: Inappropriate HDU type intercepted\n");
    return NULL;
  };
  if(phdu->pars) {
    fprintf(stderr, "new_Gpars: phdu->pars has already been allocated!\n");
  } else if(phdu->pcount > 0) {
/*
 * Allocate the new array.
 */
    phdu->pars = (Gpar *) malloc(sizeof(Gpar) * phdu->pcount);
    if(phdu->pars == NULL) {
      fprintf(stderr, "Insufficient memory to record group descriptors\n");
    } else {
      for(i=0; i<phdu->pcount; i++) { /* Initialize members */
	phdu->pars[i].ptype = NULL;
	phdu->pars[i].pscal = 0.0;
	phdu->pars[i].pzero = 0.0;
      };
    };
  };
  return phdu->pars;
}

/*.......................................................................
 * Allocate an array of primary-data-array axis descriptors in a primary
 * or image HDU.
 *
 * Input:
 *  phdu     Phdu * The PRIMARY or IMAGE extension phdu to allocate the
 *                  array for. phdu->naxis must have been set before
 *                  calling this function.
 * Output:
 *  return Imaxis * The pointer to the allocated array of descriptors
 *                  or NULL on error. This will also be assigned to
 *                  phdu->axes.
 */
static Imaxis *new_axes(Phdu *phdu)
{
  int i;
/*
 * Sanity checks.
 */
  if(phdu==NULL) {
    fprintf(stderr, "new_axes: NULL Phdu descriptor intercepted\n");
    return NULL;
  };
  if(phdu->type != F_PRIMARY && phdu->type != F_IMAGE) {
    fprintf(stderr, "new_axes: Inappropriate HDU type intercepted\n");
    return NULL;
  };
  if(phdu->axes) {
    fprintf(stderr, "new_axes: phdu->axes has already been allocated!\n");
  } else if(phdu->naxis>0) {
/*
 * Allocate the new array.
 */
    phdu->axes = (Imaxis *) malloc(sizeof(Imaxis) * phdu->naxis);
    if(phdu->axes == NULL) {
      fprintf(stderr,
	      "new_axes: Insufficient memory to record axis descriptors\n");
    } else {
      for(i=0; i<phdu->naxis; i++) { /* Initialize members */
	phdu->axes[i].ctype = NULL;
	phdu->axes[i].crpix = 0.0;
	phdu->axes[i].crval = 0.0;
	phdu->axes[i].cdelt = 0.0;
	phdu->axes[i].crota = 0.0;
      };
    };
  };
  return phdu->axes;
}

/*.......................................................................
 * Set FITS file description keywords in a primary HDU.
 * NB. All strings given here are optional. Unrequired keyword strings
 * should be sent as NULL. All strings are copied.
 *
 * Input:
 *  hdu        Hdu *  The primary header HDU descriptor.
 *  origin    char *  The origin of the FITS file.
 *  date_obs  char *  The observation date in the form "DD/MM/YY"
 *  telescop  char *  The name of the telescope used.
 *  instrume  char *  The instrument used on the given telescope.
 *  observer  char *  The name of the observer.
 *  object    char *  The source name observed.
 *  author    char *  The person who compiled the data in the file.
 *  referenc  char *  A reference to any published representation of the
 *                    data in this file.
 *  equinox double    The equinox in years for the celestial coordinate
 *                    system used in representing the data here-in.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
int setprim(Hdu *hdu, char *origin, char *date_obs, char *telescop,
	    char *instrume, char *observer, char *object, char *author,
	    char *referenc, double equinox)
{
  Phdu *phdu = (Phdu *) hdu;
/*
 * Sanity check the arguments.
 */
  if(hdu==NULL) {
    fprintf(stderr, "setprim: Intercepted NULL argument HDU descriptor\n");
    return 1;
  };
  if(hdu->complete) {
    fprintf(stderr,
	    "setprim: Illegal attempt to change established HDU structure\n");
    return 1;
  };
  if(hdu->type != F_PRIMARY) {
    fprintf(stderr, "setprim: Inappropriate HDU type intercepted\n");
    return 1;
  };
/*
 * Assign the string keyvalues taking care to avoid reassigning
 * existing values.
 */
  if(!phdu->origin && origin)
    phdu->origin = fitsstr(origin);
  if(!phdu->date_obs && date_obs)
    phdu->date_obs = fitsstr(date_obs);
  if(!phdu->telescop && telescop)
    phdu->telescop = fitsstr(telescop);
  if(!phdu->instrume && instrume)
    phdu->instrume = fitsstr(instrume);
  if(!phdu->observer && observer)
    phdu->observer = fitsstr(observer);
  if(!phdu->object && object)
    phdu->object = fitsstr(object);
  if(!phdu->author && author)
    phdu->author = fitsstr(author);
  if(!phdu->referenc && referenc)
    phdu->referenc = fitsstr(referenc);
/*
 * Other values.
 */
  if(equinox > 0.0)
    phdu->equinox = equinox;
  return 0;
}

/*.......................................................................
 * Write the header lines for the derived parts of a primary or image
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
static ADDFN(add_phdu)
{
  Phdu *phdu = (Phdu *) hdu;
  Imaxis *faxis;
  Gpar *gpar;
  int i;
/*
 * Allways signal the possible presence of extensions to the primary
 * HDU.
 */
  if(hdu->type==F_PRIMARY &&
     wlogkey(fits, hdu, "EXTEND", 0, 'T', "Extensions may be present"))
    return 1;
/*
 * Is this a random-groups primary header?
 */
  if(hdu->groups && wlogkey(fits, hdu, "GROUPS", 0, 'T', "Random-groups HDU"))
    return 1;
  if(wintkey(fits, hdu, "PCOUNT", 0, hdu->pcount, "Parameter count") ||
     wintkey(fits, hdu, "GCOUNT", 0, hdu->gcount, "Group count"))
    return 1;
/*
 * Write extension name and version keywords where given and relevant.
 */
  if(w_extkeys(fits, hdu))
    return 1;
/*
 * Write the axis description keywords where given.
 */
  if(phdu->axes) {
    for(i=0; i<phdu->naxis; i++) {
      faxis = &phdu->axes[i];
/*
 * Write the axis name.
 */
      if(faxis->ctype &&
	 wstrkey(fits, hdu, "CTYPE", i+1, faxis->ctype, "Axis name"))
	return 1;
/*
 * If the index information makes sense, write it. Note that new_axes()
 * initializes cdelt to the non-sensical value of 0.0.
 */
      if(faxis->cdelt != 0.0) {
	if(wfltkey(fits, hdu, "CRPIX", i+1, faxis->crpix, "Reference pixel") ||
	   wfltkey(fits, hdu, "CRVAL", i+1, faxis->crval, "Reference value") ||
	   wfltkey(fits, hdu, "CDELT", i+1, faxis->cdelt, "Pixel increment") ||
	   wfltkey(fits, hdu, "CROTA", i+1, faxis->crota, "Axis rotation"))
	  return 1;
      };
    };
  };
/*
 * Write the random-parameter description keywords where given.
 */
  if(phdu->pars) {
    for(i=0; i<phdu->pcount; i++) {
      gpar = &phdu->pars[i];
      if(gpar->ptype) {
	if(wstrkey(fits, hdu, "PTYPE", i+1, gpar->ptype, "Parameter name") ||
	   wfltkey(fits, hdu, "PSCAL", i+1, gpar->pscal, "Parameter scale") ||
	   wfltkey(fits, hdu, "PZERO", i+1, gpar->pzero, "Parameter offset"))
	  return 1;
      };
    };
  };
/*
 * Write the optional reserved keywords of the primary header where given.
 */
  if(hdu->type == F_PRIMARY) {
    if(phdu->origin && wstrkey(fits, hdu, "ORIGIN", 0, phdu->origin,
			       "Origin of data"))
      return 1;
    if(phdu->date_obs && wstrkey(fits, hdu, "DATE-OBS", 0, phdu->date_obs,
				 "Observation date"))
      return 1;
    if(phdu->telescop && wstrkey(fits, hdu, "TELESCOP", 0, phdu->telescop,
				 "Telescope used"))
      return 1;
    if(phdu->instrume && wstrkey(fits, hdu, "INSTRUME", 0, phdu->instrume,
				 "Instrument used"))
      return 1;
    if(phdu->observer && wstrkey(fits, hdu, "OBSERVER", 0, phdu->observer,
				 "Observers name"))
      return 1;
    if(phdu->object && wstrkey(fits, hdu, "OBJECT", 0, phdu->object,
				 "Name of observed source"))
      return 1;
    if(phdu->author && wstrkey(fits, hdu, "AUTHOR", 0, phdu->author,
				 "The author of this data"))
      return 1;
    if(phdu->referenc && wstrkey(fits, hdu, "REFERENC", 0, phdu->referenc,
				 "Published source of data"))
      return 1;
    if(phdu->equinox && wfltkey(fits, hdu, fits->aips ? "EPOCH":"EQUINOX",
				0, phdu->equinox, "Equinox of coordinates"))
      return 1;    
  };
/*
 * Write the optional reserved keywords that are common to both the
 * primary header and the IMAGE extension.
 */
  if(wfltkey(fits, hdu, "BSCALE", 0, phdu->bscale, "Scale factor of array") ||
     wfltkey(fits, hdu, "BZERO",  0, phdu->bzero,  "Zero offset of array"))
    return 1;
  if(phdu->bunit && wstrkey(fits, hdu, "BUNIT", 0, phdu->bunit,
				 "Unit of measurement"))
    return 1;
/*
 * Was a magic blanking value provided?
 */
  if(phdu->blank!=NONULL &&
     wintkey(fits, hdu, "BLANK", 0, phdu->blank, "NULL value in array"))
    return 1;
/*
 * Were min and max data values provided?
 */
  if(phdu->datamax != phdu->datamin || phdu->datamax != 0.0) {
    if(wfltkey(fits, hdu, "DATAMIN", 0, phdu->datamin, "Min data value") ||
       wfltkey(fits, hdu, "DATAMAX", 0, phdu->datamax, "Max data value"))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Copy the derived parts of a Phdu descriptor into another Phdu
 * descriptor. This function should be called only via copy_Hdu().
 *
 * Input:
 *  hdu    Hdu * The original hdu to be copied.
 * Output:
 *  return Hdu * The new copied version of 'hdu', or NULL on error.
 */
static COPYFN(cop_phdu)
{
  Phdu *old = (Phdu *) hdu;
  Hdu *new;      /* The new HDU descriptor */
  Imaxis *axis;  /* An axis descriptor from the original HDU descriptor */
  Gpar *par;     /* A parameter descriptor from the original HDU descriptor */
  int primary;   /* True if the HDU is of a primary HDU */
  int i;
/*
 * Primary or image HDU?
 */
  primary = old->type==F_PRIMARY;
/*
 * Create the new HDU.
 */
  if(primary) {
    new = new_primary(old->bitpix, old->naxis, old->dims, old->groups,
		      old->pcount, old->gcount);
  } else {
    new = new_image(old->bitpix, old->naxis, old->dims, old->extname, 0,
		    old->extlevel);
  };
  if(new==NULL)
    return new;
/*
 * Copy primary HDU specific info.
 */
  if(primary) {
    if(setprim(new, old->origin, old->date_obs, old->telescop, old->instrume,
		 old->observer, old->object, old->author, old->referenc,
		 old->equinox))
      return del_Hdu(new);
/*
 * Copy the group-parameter descriptors.
 */
    if(old->groups && old->pars) {
      for(par=old->pars,i=1; i<=old->pcount; i++,par++) {
	if(setgroup(new, i, par->ptype, par->pscal, par->pzero))
	  return del_Hdu(new);
      };
    };
  };
/*
 * Now the generic members.
 */
  if(setimage(new, old->bscale, old->bzero, old->bunit, old->blank,
	      old->datamin, old->datamax))
    return del_Hdu(new);
/*
 * Copy the axis descriptors.
 */
  if(old->axes) {
    for(axis=old->axes,i=1; i<=old->naxis; i++,axis++) {
      if(setaxis(new, i, axis->ctype, axis->crpix, axis->crval, axis->cdelt,
		 axis->crota))
	return del_Hdu(new);
    };
  };
  return new;
}

/*.......................................................................
 * Find an image axis by name.
 *
 * Input:
 *  phdu   Phdu *  The primary or IMAGE descriptor of the HDU to be
 *                 searched.
 *  ctype  char *  The name of the axis to be looked up.
 *  fixlen  int    If > 0, then this defines the maximum number of characters
 *                 to be compared. This makes it possible to search by
 *                 prefixes.
 *  start   int    The first axis to check (1-relative axis index).
 * Output:
 *  return  int    1-relative axis number 1->naxis. If not found, 0
 *                 is returned.
 */
int find_axis(Phdu *phdu, const char *ctype, int fixlen, int start)
{
  int axis;
  Imaxis *ax;
/*
 * Sanity checks.
 */
  if(phdu==NULL || (phdu->type != F_PRIMARY && phdu->type != F_IMAGE)) {
    fprintf(stderr, "find_axis: Inappropriate HDU type intercepted\n");
    return 0;
  };
/*
 * Conver the 1-relative start index to a 0-relative index.
 */
  start--;
/*
 * Seach for an axis of the name given.
 */
  if(start < phdu->naxis && start >= 0 && phdu->axes) {
    ax = &phdu->axes[start];
    for(axis=start; axis<phdu->naxis; axis++,ax++) {
      if(ax->ctype != NULL && matchstr(ax->ctype, ctype, fixlen))
	return axis + 1;
    };
  };
  return 0;
}


/*.......................................................................
 * Find an primary HDU ranodm-groups parameter by name.
 *
 * Input:
 *  phdu   Phdu *  The primary or IMAGE descriptor of the HDU to be
 *                 searched.
 *  ptype  char *  The name of the parameter to be looked up.
 *  fixlen  int    If > 0, then this defines the maximum number of characters
 *                 to be compared. This makes it possible to search by
 *                 prefixes.
 *  start   int    The first random parameter to check (1-relative index).
 * Output:
 *  return  int    1-relative parameter number 1->pcount. If not found, 0
 *                 is returned.
 */
int find_gpar(Phdu *phdu, const char *ptype, int fixlen, int start)
{
  int ipar;
  Gpar *par;
/*
 * Sanity checks.
 */
  if(phdu==NULL || phdu->type != F_PRIMARY) {
    fprintf(stderr, "find_gpar: Inappropriate HDU type intercepted\n");
    return 0;
  };
/*
 * Conver the 1-relative start index to a 0-relative index.
 */
  start--;
/*
 * Seach for an parameter of the name given.
 */
  if(start < phdu->pcount && start >= 0 && phdu->pars) {
    par = &phdu->pars[start];
    for(ipar=start; ipar<phdu->pcount; ipar++,par++) {
      if(par->ptype!=NULL && matchstr(par->ptype, ptype, fixlen))
	return ipar + 1;
    };
  };
  return 0;
}

/*.......................................................................
 * Return the name of a given random-group parameter.
 *
 * Input:
 *  phdu   Phdu *  The descriptor of the random-group primary HDU.
 *  ipar    int    The 1-relative (1->pcount) index of the parameter.
 * Output:
 *  return char *  The parameter name, or NULL on error.
 */
char *gpar_name(Phdu *phdu, int ipar)
{
  if(phdu==NULL || phdu->type!=F_PRIMARY) {
    fprintf(stderr, "gpar_name: Invalid Phdu descriptor intercepted\n");
    return NULL;
  };
  if(ipar < 1 || ipar > phdu->pcount) {
    fprintf(stderr, "gpar_name: Out-of range parameter index\n");
    return NULL;
  };
  return phdu->pars ? phdu->pars[ipar-1].ptype : NULL;
}

/*.......................................................................
 * Return the descriptor of a given random-group parameter.
 *
 * Input:
 *  phdu   Phdu *  The descriptor of the random-group primary HDU.
 *  ipar    int    The 1-relative (1->pcount) index of the parameter.
 * Output:
 *  return Gpar *  The parameter descriptor, or NULL on error.
 */
Gpar *get_gpar(Phdu *phdu, int ipar)
{
  if(phdu==NULL || phdu->type!=F_PRIMARY) {
    fprintf(stderr, "get_gpar: Invalid Phdu descriptor intercepted\n");
    return NULL;
  };
  if(ipar < 1 || ipar > phdu->pcount) {
    fprintf(stderr, "get_gpar: Out-of range parameter index\n");
    return NULL;
  };
  return phdu->pars ? &phdu->pars[ipar-1] : NULL;
}

/*.......................................................................
 * Return the axis'th axis descriptor of an IMAGE or PRIMARY HDU.
 *
 * Input:
 *  phdu      Phdu *  The descriptor of the PRIMARY or IMAGE HDU.
 *  axis       int    The (1-relative) axis (1->naxis) to return details
 *                    about.
 * Output:
 *  return  Imaxis *  The pointer to the appropriate axis descriptor, or
 *                    NULL on error.
 */
Imaxis *get_axis(Phdu *phdu, int axis)
{
/*
 * Sanity checks.
 */
  if(phdu==NULL || (phdu->type != F_PRIMARY && phdu->type != F_IMAGE)) {
    fprintf(stderr, "get_axis: Inappropriate HDU type intercepted\n");
    return NULL;
  };
/*
 * Check axis bounds.
 */
  if(axis < 1 || axis > phdu->naxis) {
    fprintf(stderr, "get_axis: 'axis' index out of bounds\n");
    return NULL;
  };
/*
 * Return the axis descriptor.
 */
  return (phdu->axes) ? &phdu->axes[axis-1] : NULL;
}

/*.......................................................................
 * Return the name of a given axis.
 *
 * Input:
 *  phdu   Phdu *  The descriptor of a primary HDU.
 *  axis     int    The 1-relative (1->naxis) index of the axis.
 * Output:
 *  return char *  The axis name, or NULL on error.
 */
char *axis_name(Phdu *phdu, int axis)
{
  if(phdu==NULL || phdu->type!=F_PRIMARY) {
    fprintf(stderr, "axis_name: Invalid Phdu descriptor intercepted\n");
    return NULL;
  };
  if(axis < 1 || axis > phdu->naxis) {
    fprintf(stderr, "axis_name: Out-of range axis index\n");
    return NULL;
  };
  return phdu->axes ? phdu->axes[axis-1].ctype : NULL;
}

/*.......................................................................
 * Complete the data section of a primary HDU.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file in which the HDU resides.
 *  hdu      Hdu *  The descriptor of the HDU to be completed.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static ENDFN(end_phdu)
{
  int waserr=0; /* Error status */
/*
 * If this is a random-groups HDU, then re-write the GCOUNT
 * header keyword with the now known number of groups.
 */
  if(hdu->groups) {
    int saveline; /* The current header line number */
    Fitkey key;   /* The descriptor of the existing GCOUNT keyword */
/*
 * Locate the GCOUNT header line.
 */
    saveline = new_hline(hdu, 0);
/*
 * Find and read the existing GCOUNT header line.
 */
    waserr = get_key(fits, hdu, "GCOUNT", DAT_INT, EOH_SEEK, &key);
/*
 * Over-write it if it needs to be changed.
 */
    if(!waserr && KEYINT(key) != hdu->gcount) {
      new_hline(hdu, what_hline(hdu)-1);
      waserr=wintkey(fits, hdu, "GCOUNT", 0, hdu->gcount, "Number of groups.");
    };
/*
 * Restore the header-line pointer.
 */
    new_hline(hdu, saveline);
  };
  return waserr!=0;
}
