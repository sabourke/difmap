#include <stdio.h>
#include <stdlib.h>

#include "logio.h"
#include "obs.h"
#include "mapmem.h"
#include "model.h"
#include "vlbutil.h"
#include "vlbconst.h"
#include "libfits.h"
#include "wmap.h"

static int primhdu(Fits *fits, Observation *ob, MapBeam *mb, int domap);
static int primdata(Fits *fits, Phdu *phdu, Observation *ob, MapBeam *mb,
		   int domap);
static int cctable(Fits *fits, Model *mod);
static int fitserr(Fits *fits);
static int hduerr(Hdu *hdu);
static int enddata(float *rowbuf, int iret);
static char *projkeyword(Observation *ob, char *key);

/*.......................................................................
 * Write the map or beam image from a MapBeam container to a FITS file.
 * Also write any associated model in a binary CC extension table.
 *
 * Input:
 *  ob  Observation *   The descriptor of the observation from which the
 *                      map/beam was derived.
 *  mb      MapBeam *   The container of the map or beam to be written.
 *  domap       int     If true save the map in the FITS file. If false
 *                      (0) save the dirty beam in the FITS file.
 *  fname      char *   The name to give the new FITS file.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
int w_MapBeam(Observation *ob, MapBeam *mb, int domap, const char *fname)
{
  Fits *fits;   /* The FITS-file descriptor */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "w_MapBeam"))
    return 1;
  if(mb==NULL) {
    lprintf(stderr, "w_MapBeam: NULL MapBeam container intercepted.\n");
    return 1;
  };
  if(fname==NULL) {
    lprintf(stderr, "w_MapBeam: NULL file name string intercepted.\n");
    return 1;
  };
/*
 * Attempt to open the new FITS file.
 */
  fits = new_Fits(fname, 0, 0, 0, 1);
  if(fits==NULL)
    return 1;
/*
 * Progress report.
 */
  lprintf(stdout, "Writing %s %s to FITS file: %s\n",
	  (domap && mb->ncmp>0) ? "clean":"dirty", domap ? "map":"beam", fname);
/*
 * Write the primary HDU containing the map or beam.
 */
  if(primhdu(fits, ob, mb, domap))
    return fitserr(fits);
/*
 * Write the CC table.
 */
  if(domap && cctable(fits, ob->model))
    return fitserr(fits);
/*
 * Close the completed FITS file and release its descriptor.
 */
  fits = del_Fits(fits);
  return 0;
}

/*.......................................................................
 * Private error cleanup funtion. This deletes the corrupt FITS descriptor.
 * and returns 1.
 */
static int fitserr(Fits *fits)
{
  fits=del_Fits(fits);
  return 1;
}

/*.......................................................................
 * Private error cleanup funtion. This deletes a given corrupt HDU
 * descriptor and returns 1. Do not use with an HDU that has been added
 * to a FITS file.
 */
static int hduerr(Hdu *hdu)
{
  hdu=del_Hdu(hdu);
  return 1;
}

/*.......................................................................
 * Create the primary HDU of the fits file in which the given map or beam
 * is to be written, and construct its header.
 *
 * Input:
 *  fits       Fits *   The FITS-file descriptor.
 *  ob  Observation *   The descriptor of the observation from which the
 *                      map/beam was derived.
 *  mb      MapBeam *   The container of the map or beam to be written.
 *  domap       int     If true save the map in the FITS file. If false
 *                      (0) save the dirty beam in the FITS file.
 * Output:
 *  return      int     0 - OK.
 */
static int primhdu(Fits *fits, Observation *ob, MapBeam *mb, int domap)
{
  enum {NAXIS=4};  /* The number of image axes */
  int dims[NAXIS]; /* This will record the dimensions of each axis */
  Hdu *hdu;        /* The generic descriptor of the primary HDU */
  int ixa,ixb;     /* Start and end X indexes in input image 0-relative */
  int iya,iyb;     /* Start and end Y indexes in input image 0-relative */
  float datamin;   /* The min value in the map or beam */
  float datamax;   /* The max value in the map or beam */
  double ra,dec;   /* Right Ascension and declination of map center (radians) */
  char history[81];/* Buffer to compose history lines in */
  int i;
/*
 * Maps are half the width of beams - determine the start and
 * end indexes along each axis. (0 - relative indexing C-style).
 */
  if(domap) {
    ixa = mb->nx/4;
    ixb = 3*ixa - 1;
    iya = mb->ny/4;
    iyb = 3*iya - 1;
  } else {
    ixa = iya = 0;
    ixb = mb->nx-1;
    iyb = mb->ny-1;
  };
/*
 * Record the dimension of each axis.
 */
  dims[0] = ixb-ixa+1;  /* RA---SIN */
  dims[1] = iyb-iya+1;  /* DEC--SIN */
  dims[2] = 1;          /* FREQ */
  dims[3] = 1;          /* STOKES */
/*
 * Create the primary HDU descriptor.
 */
  hdu = new_primary(B_FLOAT, NAXIS, dims, 0, 0, 1);
  if(hdu==NULL)
    return 1;
/*
 * Determine the RA and DEC of the map center.
 */
  ra = lmtora(ob->source.ra, ob->source.dec, -ob->geom.east, -ob->geom.north,
	      ob->proj);
  dec = lmtodec(ob->source.ra, ob->source.dec, -ob->geom.east, -ob->geom.north,
		ob->proj);
/*
 * Define the image axes.
 */
  if(setaxis(hdu, 1, projkeyword(ob, "RA"), dims[0]/2, ra * rtod,
	     -mb->xinc * rtod, 0.0) ||
     setaxis(hdu, 2, projkeyword(ob, "DEC"), dims[1]/2+1, dec * rtod,
	     mb->yinc * rtod, 0.0) ||
     setaxis(hdu, 3, "FREQ", 1, getfreq(ob,-1), getbw(ob,-1), 0.0) ||
     setaxis(hdu, 4, "STOKES", 1, (double) ob->stream.pol.type, 1.0, 0.0))
    return hduerr(hdu);
/*
 * Describe the origin of the data.
 */
  if(setprim(hdu, ob->misc.origin, ob->misc.date_obs, ob->misc.telescop,
	     ob->misc.instrume, ob->misc.observer, ob->source.name,
	     NULL, NULL, ob->misc.equinox))
    return hduerr(hdu);
/*
 * Get the minimum and maximum values of the data in the specified image
 * array.
 */
  imran((domap)?mb->map:mb->beam, mb->nx, mb->ny, ixa, ixb, iya, iyb,
	&datamin, &datamax);
/*
 * Describe the image data chracteristics.
 */
  if(setimage(hdu, 1.0, 0.0, domap?"JY/BEAM":"/BEAM", NONULL, datamin, datamax))
    return hduerr(hdu);
/*
 * Add the HDU to the FITS file.
 */
  if(add_Hdu(fits, hdu))
    return hduerr(hdu);
/*
 * Write keywords to record the beam size and number of restored components
 * if this is a restored map.
 */
  if(domap && mb->ncmp) {
    if(wfltkey(fits, hdu, "BMAJ", 0, mb->bmaj*rtod,
	       "Clean beam major axis diameter (degrees).") ||
       wfltkey(fits, hdu, "BMIN", 0, mb->bmin*rtod,
	       "Clean beam minor axis diameter (degrees).") ||
       wfltkey(fits, hdu, "BPA", 0, mb->bpa*rtod,
	       "Clean beam position angle (degrees).") ||
       wintkey(fits, hdu, "NITER", 0, mb->ncmp,
	       "Number of model components."))
      return 1;
  };
/*
 * Record the observing center if known.
 */
  if(ob->source.have_obs) {
    if(wfltkey(fits, hdu, "OBSRA", 0, ob->source.obsra * rtod,
	       "Antenna pointing RA") ||
       wfltkey(fits, hdu, "OBSDEC", 0, ob->source.obsdec * rtod,
	       "Antenna pointing Dec"))
      return 1;
  };
/*
 * Record the estimated map noise.
 */
  if(wfltkey(fits, hdu, "NOISE", 0, mb->noise,
	     "Theoretical RMS noise estimate"))
    return 1;
/*
 * Write history lines to the FITS file.
 * First rewind the history scratch file.
 */
  rec_rewind(ob->his);
  for(i=0; i<ob->nhist; i++) {
/*
 * Read the next line of history from the history scratch file.
 */
    if(rec_read(ob->his, 80, 1, history) < 0)
      return 1;
/*
 * Terminate the record.
 */
    history[80] = '\0';
/*
 * Append it to the FITS header.
 */
    if(wcomkey(fits, hdu, "HISTORY", 0, history, NULL))
      return 1;
  };
/*
 * Add an AIPS IMCLASS history line to tell AIPS what extension to give
 * the file in its catalogue.
 */
  sprintf(history, "AIPS IMCLASS=\'%s%s\'", Stokes_name(ob->stream.pol.type),
	  domap ? (mb->ncmp ? "CLN":"MAP") : "BEAM");
  if(wcomkey(fits, hdu, "HISTORY", 0, history, NULL))
    return 1;
/*
 * Add a line of history to record what has been saved.
 */
  sprintf(history, "DIFMAP  Saved %s to fits file.",
	  domap ? (mb->ncmp ? "clean-map":"residual-map"):"dirty-beam");
  if(wcomkey(fits, hdu, "HISTORY", 0, history, NULL))
    return 1;
/*
 * End the header.
 */
  if(end_header(fits, hdu))
    return 1;
/*
 * Write the primary data array.
 */
  if(primdata(fits, (Phdu *)hdu, ob, mb, domap) || end_data(fits, hdu))
    return 1;
  return 0;
}

/*.......................................................................
 * Private function of primhdu() to write a map or beam to the data
 * section of the primary HDU created by primhdu().
 *
 * Input:
 *  fits       Fits *   The FITS-file descriptor.
 *  phdu       Phdu *   The descriptor of the primary HDU.
 *  ob  Observation *   The descriptor of the observation from which the
 *                      map/beam was derived.
 *  mb      MapBeam *   The container of the map or beam to be written.
 *  domap       int     If true save the map in the FITS file. If false
 *                      (0) save the dirty beam in the FITS file.
 * Output:
 *  return      int     0 - OK.
 */
static int primdata(Fits *fits, Phdu *phdu, Observation *ob, MapBeam *mb,
		    int domap)
{
  float *rowbuf=NULL; /* Dynamically allocated buffer for row reversal */
  int ixa,ixb;     /* Start and end X indexes in input image 0-relative */
  int iya,iyb;     /* Start and end Y indexes in input image 0-relative */
  int ix,iy;       /* Indexes of pixels in the FITS file */
  int numx;        /* Length of a row in the output file */
  int numy;        /* Length of a column in the output file */
/*
 * Maps are half the width of beams - determine the start and
 * end indexes along each axis. (0 - relative indexing C-style).
 */
  if(domap) {
    ixa = mb->nx/4;
    ixb = 3*ixa - 1;
    iya = mb->ny/4;
    iyb = 3*iya - 1;
  } else {
    ixa = iya = 0;
    ixb = mb->nx-1;
    iyb = mb->ny-1;
  };
/*
 * Determine the number of pixels along the X and Y dimensions of the sub-image.
 */
  numx = ixb - ixa + 1;
  numy = iyb - iya + 1;
/*
 * Allocate a temporary buffer in which to reverse the order of pixels in
 * a row. This is required because AIPS display routines can not
 * handle the pixels being in order of increasing RA, even if CRDELTn is
 * given the correct sign.
 */
  rowbuf = malloc(numx * sizeof(float));
  if(rowbuf == NULL) {
    lprintf(stderr, "primdata: Insufficient memory for row buffer.\n");
    return 1;
  };
/*
 * Write a row of the image at a time.
 */
  for(iy=0; iy<numy; iy++) {
    float *rowptr = (domap ? mb->map:mb->beam) + mb->nx * (iy+iya);
/*
 * Reverse the new image row while copying it into rowbuf[].
 */
    for(ix=0; ix<numx; ix++)
      rowbuf[ix] = rowptr[ixb-ix];
/*
 * Write the row to the FITS file.
 */
    if(wimage(fits, phdu, 0L, (long)numx*iy, numx, DAT_FLT, 0, NULL, rowbuf) < numx)
      return enddata(rowbuf, 1);
  };
  return enddata(rowbuf, 0);
}

/*.......................................................................
 * Private cleanup function of primdata() used to free the row buffer
 * and return the given primdata() return code.
 *
 * Input:
 *  rowbuf    float *  The buffer to be freed.
 *  iret        int    The desired return code.
 * Output:
 *  return      int    iret.
 */
static int enddata(float *rowbuf, int iret)
{
  if(rowbuf)
    free(rowbuf);
  return iret;
}

/*.......................................................................
 * Write an AIPS CC table for the given model.
 *
 * Input:
 *  fits       Fits *   The descriptor of the opened FITS file.
 *  mod       Model *   The model to be written to a CC table. A CC
 *                      table will only be written only if mod!=NULL
 *                      and mod->ncmp>0.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
static int cctable(Fits *fits, Model *mod)
{
  Hdu *hdu;    /* The generic HDU descriptor of the table */
  Thdu *thdu;  /* The generic table HDU descriptor of the table */
  Modcmp *cmp; /* A model component */
  int irow;    /* The number of the row being written */
  int i;
/*
 * Define the 7 CC table columns.
 */
  static struct {char *tform; char *ttype; char *tunit;} cctab[]= {
    {"1E", "FLUX",     "JY"},
    {"1E", "DELTAX",   "DEGREES"},
    {"1E", "DELTAY",   "DEGREES"},
    {"1E", "MAJOR AX", "DEGREES"},
    {"1E", "MINOR AX", "DEGREES"},
    {"1E", "POSANGLE", "DEGREES"},
    {"1E", "TYPE OBJ", "CODE"},
  };
  static const int ncol = sizeof(cctab)/sizeof(cctab[0]);
/*
 * Is there a model to be written.
 */
  if(mod && mod->ncmp>0) {
    int ncmp=0;  /* The actual number of representable components */
/*
 * Count the number of model components that can be represented in an
 * AIPS CC table.
 */
    for(cmp=mod->head; cmp; cmp=cmp->next) {
      switch(cmp->type) {
      case M_DELT: case M_GAUS:
	ncmp++;
	break;
      default:
	break;
      };
    };
/*
 * Report un-representable components.
 */
    if(ncmp < mod->ncmp) {
      lprintf(stderr, "cctable: Warning: Only gaussian and delta function components can be\n");
      lprintf(stderr, " re-represented in AIPS CC tables.\n");
      lprintf(stderr, " For this reason %d components have been omitted.\n",
	      mod->ncmp-ncmp);
    };
/*
 * Are there any components left to be written?
 */
    if(ncmp<=0)
      return 0;
/*
 * Create the HDU descriptor of the new binary CC table.
 */
    hdu = new_bintab(ncmp, "AIPS CC", 0, 0, ncol, 0L);
    thdu = (Thdu *) hdu;
/*
 * Describe the table columns.
 */
    for(i=0;i<ncol;i++) {
      if(setbfield(hdu, i+1, 1.0, 0.0, cctab[i].tform, NONULL, cctab[i].ttype,
		   cctab[i].tunit, NULL, NULL))
	return hduerr(hdu);
    };
/*
 * Add the HDU to the FITS file.
 */
    if(add_Hdu(fits, hdu))
      return hduerr(hdu);
/*
 * Terminate the header.
 */
    if(end_header(fits, hdu))
      return 1;
/*
 * Write the model components to the table.
 */
    irow = 0;
    for(cmp=mod->head; cmp; cmp=cmp->next) {
      int cctype = -1;
/*
 * Check that the model component can be represented in a CC table and
 * get its AIPS CC type code.
 */
      switch(cmp->type) {
      case M_DELT:
	cctype = 0;
	break;
      case M_GAUS:
	cctype = 1;
	break;
      default:
	break;
      };
/*
 * Write representable components to the file, and leave the columns
 * allocated for un-representable components as zeroes (zero flux delta
 * components).
 */
      if(cctype>=0) {
	float colval;   /* The value to be written to a given column */
/*
 * New row.
 */
	irow++;
/*
 * FLUX.
 */
	colval = cmp->flux;
	if(wcolumn(fits, thdu, 1, irow, DAT_FLT, 0, NULL, 0, 1, &colval)<1)
	  return 1;
/*
 * DELTAX.
 */
	colval = cmp->x * rtod;
	if(wcolumn(fits, thdu, 2, irow, DAT_FLT, 0, NULL, 0, 1, &colval)<1)
	  return 1;
/*
 * DELTAY.
 */
	colval = cmp->y * rtod;
	if(wcolumn(fits, thdu, 3, irow, DAT_FLT, 0, NULL, 0, 1, &colval)<1)
	  return 1;
/*
 * MAJOR AX.
 */
	colval = cmp->major * rtod;
	if(wcolumn(fits, thdu, 4, irow, DAT_FLT, 0, NULL, 0, 1, &colval)<1)
	  return 1;
/*
 * MINOR AX.
 */
	colval = cmp->ratio * cmp->major * rtod;
	if(wcolumn(fits, thdu, 5, irow, DAT_FLT, 0, NULL, 0, 1, &colval)<1)
	  return 1;
/*
 * POSANGLE.
 */
	colval = cmp->phi * rtod;
	if(wcolumn(fits, thdu, 6, irow, DAT_FLT, 0, NULL, 0, 1, &colval)<1)
	  return 1;
/*
 * TYPE OBJ.
 */
	colval = cctype;
	if(wcolumn(fits, thdu, 7, irow, DAT_FLT, 0, NULL, 0, 1, &colval)<1)
	  return 1;
      };
    };
/*
 * Close the data-section of the table.
 */
    if(end_data(fits, hdu))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Combine a name of <= 4 characters with a projection code, to produce
 * a header keyword name. Return the name as an 8 character '\0' string
 * padded with hyphens.
 *
 * Input:
 *  ob    Observation *   The observation from which to determine the
 *                        projection.
 *  key          char *   The keyword name (eg. "RA" or "DEC").
 * Output:
 *  return       char *   The output keyword name.
 */
static char *projkeyword(Observation *ob, char *key)
{
  static char name[9];  /* The output string */
/*
 * Get the projection name.
 */
  char *proj = Proj_name(ob->proj);
  sprintf(name, "%-4.4s%4.4s", key, proj);
/*
 * Convert interword spaces to hyphens.
 */
  {
    char *cptr;
    for(cptr=name; *cptr; cptr++) {
      if(*cptr == ' ')
	*cptr = '-';
    };
  };
  return name;
}

