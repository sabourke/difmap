#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "logio.h"
#include "libfits.h"
#include "obs.h"
#include "vlbconst.h"

/*
 * The primary HDU axes will be arranged in the following order:
 *   COMPLEX, STOKES, FREQ, IF, RA, DEC.
 * The random group parameters will be arranged in the following order:
 *   UU, VV, WW, BASELINE, DATE, DATE
 */

enum {NAXIS=7};  /* There will always be 7 axes */

static int wrterr(Fits *fits);

static int primhdu(Observation *ob, Fits *fits, int doshift);
static int primdata(Observation *ob, Fits *fits, Phdu *phdu, int doshift);
static int prim_err(double *gpar, double *data);
static int hduerr(Hdu *hdu);

static int fqtable(Observation *ob, Fits *fits);
static int fqdata(Observation *ob, Fits *fits, Thdu *thdu);

static int antable(Fits *fits, Subarray *sub, int extver);
static int ascan(Fits *fits, Subarray *sub, int extver);
static int ascandata(Fits *fits, Thdu *thdu, Subarray *sub);
static int binan(Fits *fits, Subarray *sub, int extver);
static int binandata(Fits *fits, Thdu *thdu, Subarray *sub);
static char *uvwname(Observation *ob, char *uvw);
static int wrt_p_refant(Fits *fits, Hdu *hdu, Subarray *sub);

/*.......................................................................
 * Write a UV FITS file from the contents of an Observation structure.
 *
 * Input:
 *  ob   Observation *  The descriptor of the observation to be written.
 *  name  const char *  The name of the output UV FITS file.
 *  doshift      int    If true, apply the shift in ob->geom to the data,
 *                      and change the recorded RA and Dec accordingly.
 * Output:
 *  return       int    0 - OK.
 *                      1 - Error.
 */
int uvf_write(Observation *ob, const char *name, int doshift)
{
  Fits *fits;   /* The descriptor of the FITS file */
  int isub;     /* The sub-array being written */
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "uvf_write"))
    return 1;
  if(name==NULL) {
    lprintf(stderr, "uvf_write: NULL file name intercepted.\n");
    return 1;
  };
/*
 * Apply cached edits.
 */
  if(ed_flush(ob))
    return 1;
/*
 * Create the new FITS file.
 */
  fits = new_Fits(name, 0, 0, 0, 1);
  if(fits==NULL)
    return 1;
/*
 * Keep user informed.
 */
  lprintf(stdout, "Writing UV FITS file: %s\n", name);
/*
 * Create the primary HDU.
 */
  if(primhdu(ob, fits, doshift))
    return wrterr(fits);
/*
 * Write an AIPS FQ table.
 */
  if(fqtable(ob, fits))
    return wrterr(fits);
/*
 * Write an AIPS AN antenna table for each sub-array.
 */
  for(isub=0; isub<ob->nsub; isub++) {
    if(antable(fits, &ob->sub[isub], isub+1))
      return wrterr(fits);
  };
/*
 * Close the completed file.
 */
  fits = del_Fits(fits);
  return 0;
}

/*.......................................................................
 * Private error cleanup function of uvf_write().
 */
static int wrterr(Fits *fits)
{
  fits = del_Fits(fits);
  return 1;
}

/*.......................................................................
 * Construct and write the header of the primary HDU.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation being written.
 *  fits       Fits *  The descriptor of the FITS file.
 *  doshift     int    If true, apply the shift in ob->geom to the data,
 *                     and change the recorded RA and Dec accordingly.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int primhdu(Observation *ob, Fits *fits, int doshift)
{
  int dims[NAXIS]; /* This will record the dimensions of each axis */
  Hdu *hdu;        /* The generic descriptor of the primary HDU */
  char history[81];/* Buffer to compose history lines in */
  int npar;        /* The number of random group parameters */
  double ra,dec;   /* The right ascension and declination of the pointing */
                   /*  center. */
  int i;
/*
 * Assign the dimension of each axis.
 */
  dims[0] = 0;         /* Random groups format feature */
  dims[1] = 3;         /* COMPLEX real,imaginary,weight */
  dims[2] = ob->npol;  /* STOKES */
  dims[3] = ob->nchan; /* FREQ */
  dims[4] = ob->nif;   /* IF */
  dims[5] = 1;         /* RA */
  dims[6] = 1;         /* DEC */
/*
 * How many random group parameters do we have?
 */
  npar = ob->have_inttim ? 7 : 6;
/*
 * Create the primary HDU.
 */
  hdu = new_primary(B_FLOAT, NAXIS, dims, 1, npar, 0);
  if(hdu==NULL)
    return 1;
/*
 * Create the COMPLEX axis.
 */
  if(setaxis(hdu, 2, "COMPLEX", 1.0, 1.0, 1.0, 0.0))
    return hduerr(hdu);
/*
 * Create the STOKES axis.
 */
  if(setaxis(hdu, 3, "STOKES", 1.0, (double)ob->pols[0],
	     (ob->npol<=1 ? 1.0 : ob->pols[1]-ob->pols[0]), 0.0))
    return hduerr(hdu);
/*
 * Create the FREQ axis.
 */
  if(setaxis(hdu, 4, "FREQ", 1.0, ob->ifs[0].freq, ob->ifs[0].df, 0.0))
    return hduerr(hdu);
/*
 * Create the IF axis.
 */
  if(setaxis(hdu, 5, "IF", 1.0, 1.0, 1.0, 0.0))
    return hduerr(hdu);
/*
 * Create the RA axis.
 */
  if(doshift)
    ra = lmtora(ob->source.ra, ob->source.dec, -ob->geom.east,
		-ob->geom.north, ob->proj);
  else
    ra = ob->source.ra;
  if(setaxis(hdu, 6, "RA", 1.0, ra * rtod, 1.0, 0.0))
    return hduerr(hdu);
/*
 * Create the DEC axis.
 */
  if(doshift)
    dec = lmtodec(ob->source.ra, ob->source.dec, -ob->geom.east,
		-ob->geom.north, ob->proj);
  else
    dec = ob->source.dec;
  if(setaxis(hdu, 7, "DEC", 1.0, dec * rtod, 1.0, 0.0))
    return hduerr(hdu);
/*
 * Now initialize the random group parameters.
 *
 * Create the UU, VV and WW random group parameters.
 */
  if(setgroup(hdu, 1, uvwname(ob, "UU"), 1.0, 0.0) ||
     setgroup(hdu, 2, uvwname(ob, "VV"), 1.0, 0.0) ||
     setgroup(hdu, 3, uvwname(ob, "WW"), 1.0, 0.0))
    return hduerr(hdu);
/*
 * Create the BASELINE random group parameter.
 */
  if(setgroup(hdu, 4, "BASELINE", 1.0, 0.0))
    return hduerr(hdu);
/*
 * Create the two DATE random group parameters.
 * Split the date offset into integral and fractional parts to
 * preserve precision.
 */
  {
    double dfrc;   /* Fractional part of date */
    double dint;   /* The integral part of date (number of days) */
/*
 * Split the reference date into integral and fractional parts.
 */
    dfrc = modf(ob->date.utc_ref, &dint);
/*
 * Convert to absolute Julian date by adding 2400000.5.
 */
    dint += 2400000.0;
    dfrc += 0.5;
    if(setgroup(hdu, 5, "DATE", 1.0, dint) ||
       setgroup(hdu, 6, "DATE", 1.0, dfrc))
      return hduerr(hdu);
  };
/*
 * If available, arrange to write integration times.
 */
  if(ob->have_inttim && setgroup(hdu, 7, "INTTIM", 1.0, 0.0))
    return hduerr(hdu);
/*
 * Fill in miscellaneous details of the Observation.
 */
  if(setprim(hdu, ob->misc.origin, ob->misc.date_obs, ob->misc.telescop,
	     ob->misc.instrume, ob->misc.observer, ob->source.name, NULL,
	     NULL, ob->misc.equinox))
    return hduerr(hdu);
/*
 * Set the scaling of the data array.
 */
  if(setimage(hdu, 1.0, 0.0, ob->misc.bunit, NONULL, 0.0, 0.0))
    return hduerr(hdu);
/*
 * Add the initialized HDU to the FITS file.
 */
  if(add_Hdu(fits, hdu))
    return hduerr(hdu);
/*
 * Write velocity info if known.
 */
  if(ob->vel.velref != 0) {
    if(wintkey(fits, hdu, "VELREF", 0, ob->vel.velref,
	       ">256 RADIO, 1 LSR 2 HEL 3 OBS") ||
       wfltkey(fits, hdu, "ALTRVAL", 0, ob->vel.altrval,
	       "Alternate Freq/vel ref value") ||
       wfltkey(fits, hdu, "ALTRPIX", 0, ob->vel.altrpix,
	       "Alternate Freq/vel ref pixel") ||
       wfltkey(fits, hdu, "RESTFREQ", 0, ob->vel.restfreq,
	       "Rest frequency"))
      return hduerr(hdu);
  };
/*
 * Write the antenna pointing-center if known.
 */
  if(ob->source.have_obs) {
    wfltkey(fits, hdu, "OBSRA", 0, ob->source.obsra * rtod,
	    "Antenna pointing RA");
    wfltkey(fits, hdu, "OBSDEC", 0, ob->source.obsdec * rtod,
	    "Antenna pointing Dec");
  };
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
 * Append a special AIPS IMCLASS history line.
 */
  if(wcomkey(fits, hdu, "HISTORY", 0, "AIPS IMCLASS=\'UVF\'", NULL))
    return 1;
/*
 * Construct and append a special AIPS WTSCAL history line.
 */
  sprintf(history, "AIPS WTSCAL=%G", ob->geom.wtscale);
  if(wcomkey(fits, hdu, "HISTORY", 0, history, NULL))
    return 1;
/*
 * Append a special AIPS SORT ORDER history line.
 */
  if(wcomkey(fits, hdu, "HISTORY", 0, "AIPS SORT ORDER=\'TB\'", NULL))
    return 1;
/*
 * End the header.
 */
  if(end_header(fits, hdu))
    return 1;
/*
 * Write the data section of the primary HDU.
 */
  if(primdata(ob, fits, (Phdu *)hdu, doshift) || end_data(fits, hdu))
    return 1;
  return 0;
}

/*.......................................................................
 * Private error cleanup function used to delete a given un-installed
 * HDU descriptor before returning the error code, 1. This function
 * must not be used with an HDU that has been added to a FITS file.
 * Use like:
 *
 *    if(error)
 *      return hduerr(hdu);
 */
static int hduerr(Hdu *hdu)
{
  hdu = del_Hdu(hdu);
  return 1;
}

/*.......................................................................
 * Write the corrected UV data into the random-groups structure of the
 * primary HDU.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation being written.
 *  fits       Fits *  The descriptor of the FITS file.
 *  phdu       Phdu *  The descriptor of the primary HDU.
 *  doshift     int    If true, shift the data by the amount recorded
 *                     in ob->geom.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int primdata(Observation *ob, Fits *fits, Phdu *phdu, int doshift)
{
  Intrec *rec;      /* The integration file-record descriptor */
  long group=0L;    /* The group being written */
  long ndata;       /* Length of one group data array (double's) */
  int ut;           /* The index of the integration-record being written */
  int base;         /* The index of the baseline being written */
  int cif;          /* The index of the IF being processed */
  int fc;           /* The index of the spectral-line channel being processed */
  int pol;          /* The index of the stokes parameter being processed */
  int npar;         /* The number of random group parameters */
  double *data=NULL;/* Buffer used for composing output data in */
  double *gpar=NULL;/* Buffer used for composing random-group parameters in */
/*
 * Make sure that all defered edits have been applied to the
 * uvdata scratch file.
 */
  if(ed_flush(ob))
    return 1;
/*
 * How big is one group data-array in double's?
 */
  ndata = 3L * ob->npol * ob->nchan * ob->nif;
/*
 * How many random group parameters do we have?
 */
  npar = ob->have_inttim ? 7 : 6;
/*
 * Allocate two output buffers - one to compose random-group parameters in and
 * a second to compose the complex visibilities in.
 */
  gpar = malloc(sizeof(double) * npar);
  data = malloc(sizeof(double) * ndata);
  if(gpar==NULL || data==NULL) {
    lprintf(stderr, "primdata: Insufficient memory to write UV data.\n");
    return prim_err(gpar, data);
  };
/*
 * Initialize to read all data in the uvdata scratch file.
 */
  {
    Dpage *dp = ob->dp;
    if(dp_crange(dp, 0, ob->nchan-1) ||
       dp_irange(dp, 0, ob->nif-1)   ||
       dp_brange(dp, 0, ob->nbmax-1) ||
       dp_srange(dp, 0, ob->npol-1))
      return prim_err(gpar, data);
  };
/*
 * Loop to write in TB order.
 */
  rec = ob->rec;
  for(ut=0; ut<ob->nrec; ut++, rec++) {
    Integration *integ = rec->integ;
    Visibility *vis = integ->vis;
    Subarray *sub = integ->sub;
    Baseline *bptr = sub->base;
    Station *tel = sub->tel;
    int isub;  /* Index of the sub-array of the current integration */
/*
 * Determine the index of the sub-array for this integration.
 */
    for(isub=0; isub<ob->nsub && sub != &ob->sub[isub]; isub++);
/*
 * Read the next integration of data from the uvdata scratch file.
 */
    if(dp_read(ob->dp, (long) ut))
      return prim_err(gpar, data);
/*
 * There is one group to be written per baseline.
 */
    for(base=0; base<sub->nbase; base++,vis++,bptr++) {
      int ta = bptr->tel_a;  /* First telescope on the current baseline */
      int tb = bptr->tel_b;  /* Second telescope on the current baseline */
      int ok=0;              /* 0 - if all visibilities are deleted */
      double *dptr = data;    /* Pointer into data[] */
      Dif *ifs = ob->dp->ifs; /* Pointer to tree of visibilities in buffer */
/*
 * Compose the visibility data to be written to the file in this group.
 * Extract the data from the read buffer and copy into the output buffer.
 */
      for(cif=0; cif<ob->nif; cif++,ifs++) {
	Dchan *dchan=ifs->chan;
	for(fc=0; fc<ob->nchan; fc++,dchan++) {
	  Cvis *cvis = dchan->base[base].pol;
	  for(pol=0; pol<ob->npol; pol++,cvis++) {
	    *dptr++ = cvis->re;
	    *dptr++ = cvis->im;
	    *dptr++ = cvis->wt;
	    ok = ok || cvis->wt != 0.0f;   /* Visibility not deleted? */
	  };
	};
      };
/*
 * If any usable data was accumulated in the output buffer, apply self-cal
 * and resoff corrections to them and write them to the FITS file.
 */
      if(ok) {
	dptr = data;
	for(cif=0; cif<ob->nif; cif++) {
	  Telcor *tcor = integ->icor[cif].tcor;
	  Bascor *bcor = &bptr->bcor[cif];
	  If *ifp = ob->ifs + cif;
/*
 * Combine selfcal and resoff contributions to the amp and phase corrections.
 */
	  float amp_cor = tcor[ta].amp_cor * tcor[tb].amp_cor * bcor->amp_cor;
	  float phs_cor = tcor[ta].phs_cor - tcor[tb].phs_cor + bcor->phs_cor;
/*
 * Pre-compute the cos and sin of the phase correction for use in
 * correcting the complex representation of the visibilities.
 */
	  float cosphi = cos(phs_cor);
	  float sinphi = sin(phs_cor);
/*
 * Determine whether the correction for this baseline is flagged.
 */
	  int bad_cor = tcor[ta].bad || tcor[tb].bad;
/*
 * Ensure that the amplitude correction is +ve.
 */
	  amp_cor = fabs(amp_cor);
/*
 * All spectral-line channels and polarizations receive the same corrections
 * for a given IF (unless doshift is enabled).
 */
	  for(fc=0; fc<ob->nchan; fc++) {
/*
 * If any shifts are to be applied, modify the phase correction to
 * accomplish this.
 */
	    if(doshift) {
/*
 * We will need to evaluate the fourier component 2.pi.u.dx + 2.pi.v.dy,
 * where u and v have been converted from light seconds to wavelengths.
 * Compute the center frequency of the latest spectral line channel.
 */
	      float freq = ifp->freq + fc * ifp->df;
/*
 * Compute the phase shift needed.
 */
	      float phi = phs_cor + twopi * freq * (ob->geom.east  * vis->u +
						    ob->geom.north * vis->v);
	      cosphi = cos(phi);
	      sinphi = sin(phi);
	    };
/*
 * Apply corrections to each of the recorded polarizations.
 */
	    for(pol=0; pol<ob->npol; pol++) {
	      float re = dptr[0];
	      float im = dptr[1];
	      float wt = dptr[2];
/*
 * The complex correction is: (x+iy) * amp_cor * exp(i.phs_cor).
 */
	      *dptr++ = amp_cor * (re * cosphi - im * sinphi); /* Real */
	      *dptr++ = amp_cor * (re * sinphi + im * cosphi); /* Imaginary */
	      *dptr++ = (bad_cor && wt>0.0f ? -wt:wt) / (amp_cor * amp_cor);
	                /* Weight = 1/amp_err^2 */
	    };
	  };
	};
/*
 * Construct the random parameters for the new baseline.
 */
	gpar[0] = vis->u;
	gpar[1] = vis->v;
	gpar[2] = vis->w;
	gpar[3] = 256 * tel[ta].antno + tel[tb].antno + 0.01 * isub;
/*
 * Given that the MJD date of the start of the year of the observation
 * is stored in the POFF keyword for the first of two DATE random
 * parameters, the MJD date random parameters when combined, represent
 * the number of days since the start of the that year. Since we
 * require accuracy in that date to better than a second, and should
 * allow combined observations that span hundreds of days, while an
 * IEEE float only caters for a little over 6 decimal significant
 * figures, the date will be split into its fractional and integral
 * parts and stored in two random parameters.
 */
	{
	  double date;   /* The TAI MJD wrt the first integration */
	  double dfrc;   /* Fractional part of date */
	  double dint;   /* The integral part of date */
	  date = (integ->ut + sub->datutc) / daysec;
/*
 * Split the date into integral and fractional parts.
 */
	  dfrc = modf(date, &dint);
/*
 * Store the integral part in the first DATE parameter and the fractional
 * part in the second.
 */
	  gpar[4] = dint;
	  gpar[5] = dfrc;
	};
/*
 * If known, include the integration time.
 */
	if(ob->have_inttim)
	  gpar[6] = vis->dt;
/*
 * Write them to the FITS file, without removing scale and offsets,
 * since they have already been taken care of above.
 */
	if(wgroup(fits, phdu, group, 0L, (long) npar, DAT_DBL, 0, NULL, gpar)
	   != npar)
	  return prim_err(gpar, data);
/*
 * Write data to the FITS file.
 */
	if(wimage(fits, phdu, group, 0L, ndata, DAT_DBL, 1, NULL, data)!=ndata)
	  return prim_err(gpar, data);
/*
 * Pre-pare for the next un-deleted group.
 */
	group++;
      };
    };
  };
/*
 * Release the buffer arrays.
 */
  if(gpar)
    free(gpar);
  if(data)
    free(data);
  return 0;
}

/*.......................................................................
 * Private error cleanup function of primdata().
 */
static int prim_err(double *gpar, double *data)
{
/*
 * Release memory allocated to the group-parameter and data arrays.
 */
  if(gpar)
    free(gpar);
  if(data)
    free(data);
  return 1;
}

/*.......................................................................
 * Construct and write an AIPS FQ table.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation being written.
 *  fits       Fits *  The descriptor of the FITS file.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int fqtable(Observation *ob, Fits *fits)
{
  enum {NCOL=5};   /* The number of columns in the table */
  enum {IFMAG=10}; /* Only make allowances for up to 10^IFMAG IFs */
  char tform[IFMAG+2]; /* Room for the IF count and format letter */
  Hdu *hdu;        /* The generic descriptor of the table HDU */
/*
 * Create the HDU descriptor of a binary table.
 */
  hdu = new_bintab(1, "AIPS FQ", 1, 1, NCOL, 0L);
  if(hdu==NULL)
    return 1;
/*
 * A number of items in the table are arrays of dimension ob->nif.
 * The format of these columns is thus the value of ob->nif plus a
 * type letter. Check that the tform[] array is big enough to handle this.
 */
  if(ob->nif > pow(10.0, (double) IFMAG)) {
    lprintf(stderr, "fqtable: There are more IFs than the max expected.\n");
    return hduerr(hdu);
  };
/*
 * Describe the details of each column.
 */
  if(setbfield(hdu, 1, 1.0, 0.0, "1J", NONULL, "FRQSEL", "", NULL, NULL))
    return hduerr(hdu);
  sprintf(tform, "%dD", ob->nif);
  if(setbfield(hdu, 2, 1.0, 0.0, tform, NONULL, "IF FREQ", "HZ", NULL, NULL))
    return hduerr(hdu);
  sprintf(tform, "%dE", ob->nif);
  if(setbfield(hdu, 3, 1.0, 0.0, tform, NONULL, "CH WIDTH", "HZ", NULL, NULL))
    return hduerr(hdu);
  sprintf(tform, "%dE", ob->nif);
  if(setbfield(hdu, 4, 1.0, 0.0, tform, NONULL, "TOTAL BANDWIDTH", "HZ",NULL,NULL))
    return hduerr(hdu);
  sprintf(tform, "%dJ", ob->nif);
  if(setbfield(hdu, 5, 1.0, 0.0, tform, NONULL, "SIDEBAND", "", NULL, NULL))
    return hduerr(hdu);
/*
 * Add the initialized HDU to the FITS file.
 */
  if(add_Hdu(fits, hdu))
    return hduerr(hdu);
/*
 * Write the keyword that parameterises the number of IFs.
 */
  if(wintkey(fits, hdu, "NO_IF", 0, ob->nif, "Number of IFs"))
    return 1;
/*
 * End the header.
 */
  if(end_header(fits, hdu))
    return 1;
/*
 * Write the data section of the table.
 */
  if(fqdata(ob, fits, (Thdu *) hdu) || end_data(fits, hdu))
    return 1;
  return 0;
}

/*.......................................................................
 * Write the contents of an AIPS FQ table.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation being written.
 *  fits       Fits *  The descriptor of the FITS file.
 *  thdu       Thdu *  The descriptor of the FQ table HDU.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int fqdata(Observation *ob, Fits *fits, Thdu *thdu)
{
  If *ifptr;       /* Pointer to the descriptor of the IF being recorded */
  int cif;         /* The index of the IF being recorded */
  int itmp;        /* Used to compile integer data items in */
  double dtmp;     /* Used to compile double data items in */
/*
 * The first column contains the frequency ID associated with the row.
 * Give this the value 1.
 */
  itmp = 1;
  if(wcolumn(fits, thdu, 1, 1, DAT_INT, 0, NULL, 0, 1, &itmp)==0)
    return 1;
/*
 * Write the required members of the table row for each IF.
 */
  ifptr = &ob->ifs[0];
  for(cif=0; cif<ob->nif; cif++,ifptr++) {
/*
 * Record the IF frequency offset wrt the frequency in the first IF.
 */
    dtmp = ifptr->freq - ob->ifs[0].freq;
    if(wcolumn(fits, thdu, 2, 1, DAT_DBL, 0, NULL, cif, 1, &dtmp)==0)
      return 1;
/*
 * Record the spectral-line channel width in this IF.
 */
    dtmp = fabs(ifptr->df);
    if(wcolumn(fits, thdu, 3, 1, DAT_DBL, 0, NULL, cif, 1, &dtmp)==0)
      return 1;
/*
 * Record the total bandwidth of the IF.
 */
    dtmp = ifptr->bw;
    if(wcolumn(fits, thdu, 4, 1, DAT_DBL, 0, NULL, cif, 1, &dtmp)==0)
      return 1;
/*
 * Encode the sideband type recorded in this IF.
 */
    itmp = ifptr->df < 0.0 ? -1 : 1;
    if(wcolumn(fits, thdu, 5, 1, DAT_INT, 0, NULL, cif, 1, &itmp)==0)
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Construct and write an AIPS AN antenna table.
 *
 * Input:
 *  fits       Fits *  The descriptor of the FITS file.
 *  sub    Subarray *  The descriptor of the sub-array to be recorded.
 *  extver      int    The extension version number for the AN table.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int antable(Fits *fits, Subarray *sub, int extver)
{
/*
 * Write out the same form of AN table as was read from the original
 * FITS file for the given sub-array. The different forms contain
 * different types of info, and in particular the ASCII form contains
 * insufficient info to create the binary form.
 */
  if(sub->binan) {
    if(binan(fits, sub, extver))
      return 1;
  } else {
    if(ascan(fits, sub, extver))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Construct and write an ASCII AIPS AN antenna table.
 *
 * Input:
 *  fits       Fits *  The descriptor of the FITS file.
 *  sub    Subarray *  The descriptor of the sub-array to be recorded.
 *  extver      int    The extension version number for the AN table.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int ascan(Fits *fits, Subarray *sub, int extver)
{
  enum {NCOL=5};   /* The number of columns in the table */
  Hdu *hdu;        /* The generic descriptor of the HDU */
/*
 * Create the HDU descriptor of a binary table.
 */
  hdu = new_asctab(80, sub->nstat, "AIPS AN", extver, 1, NCOL);
  if(hdu==NULL)
    return 1;
/*
 * Describe the details of each column.
 */
  if(setafield(hdu, 1,  1, 1.0, 0.0, "I3",     NULL, "ANT NO.", "")  ||
     setafield(hdu, 2,  7, 1.0, 0.0, "A8",     NULL, "STATION", "") ||
     setafield(hdu, 3, 15, 1.0, 0.0, "D20.10", NULL, "LX", "METERS")  ||
     setafield(hdu, 4, 35, 1.0, 0.0, "D20.10", NULL, "LY", "METERS")  ||
     setafield(hdu, 5, 55, 1.0, 0.0, "D20.10", NULL, "LZ", "METERS"))
    return hduerr(hdu);
/*
 * Add the initialized HDU to the FITS file.
 */
  if(add_Hdu(fits, hdu))
    return hduerr(hdu);
/*
 * Write polarization P_REFANT and P_DIFFnn keywords if necessary.
 */
  if(wrt_p_refant(fits, hdu, sub))
    return 1;
/*
 * End the header.
 */
  if(end_header(fits, hdu))
    return 1;
/*
 * Write the data part of the ASCII AIPS AN table.
 */
  if(ascandata(fits, (Thdu *)hdu, sub) || end_data(fits, hdu))
    return 1;
  return 0;
}

/*.......................................................................
 * Write the data portion of an ASCII AIPS AN antenna table.
 *
 * Input:
 *  fits       Fits *  The descriptor of the FITS file.
 *  thdu       Thdu *  The descriptor of the table HDU.
 *  sub    Subarray *  The descriptor of the sub-array to be recorded.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int ascandata(Fits *fits, Thdu *thdu, Subarray *sub)
{
  Station *tel;  /* Pointer to the descriptor of the station being recorded */
  int irow;      /* The row in the table (telescope index + 1) */
  int itmp;      /* Used to compile integer data items in */
  double dtmp;   /* Used to compile double data items in */
/*
 * Get the pointer to the array of telescope descriptors.
 */
  tel = sub->tel;
/*
 * Each antenna is recorded in one row of the table.
 */
  for(irow=1; irow<=sub->nstat; irow++,tel++) {
/*
 * Station number.
 */
    itmp = tel->antno;
    if(wcolumn(fits, thdu, 1, irow, DAT_INT, 0, NULL, 0, 1, &itmp)==0)
      return 1;
/*
 * Station name.
 */
    itmp = strlen(tel->name);
    if(wcolumn(fits, thdu, 2, irow, DAT_CHR, 0, NULL, 0, itmp, tel->name)==0)
      return 1;
/*
 * LX.
 */
    dtmp = tel->geo.gnd.x;
    if(wcolumn(fits, thdu, 3, irow, DAT_DBL, 0, NULL, 0, 1, &dtmp)==0)
      return 1;
/*
 * LY.
 */
    dtmp = tel->geo.gnd.y;
    if(wcolumn(fits, thdu, 4, irow, DAT_DBL, 0, NULL, 0, 1, &dtmp)==0)
      return 1;
/*
 * LZ.
 */
    dtmp = tel->geo.gnd.z;
    if(wcolumn(fits, thdu, 5, irow, DAT_DBL, 0, NULL, 0, 1, &dtmp)==0)
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Construct and write a binary AIPS AN antenna table.
 *
 * Input:
 *  fits       Fits *  The descriptor of the FITS file.
 *  sub    Subarray *  The sub-array to be recorded.
 *  extver      int    The extension version number for the AN table.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int binan(Fits *fits, Subarray *sub, int extver)
{
  enum {NCOL=12};  /* The number of columns in the table */
  char tform[20];  /* Array to compile table column formats in */
  Hdu *hdu;        /* The generic descriptor of the HDU */
  Binan *an;       /* The container of binary AN table info. */
/*
 * Get the container of input binary AN-table info.
 */
  an = sub->binan;
/*
 * Create the HDU descriptor of a binary table.
 */
  hdu = new_bintab(sub->nstat, "AIPS AN", extver, 1, NCOL, 0L);
  if(hdu==NULL)
    return 1;
/*
 * Describe the details of each column.
 */
  if(setbfield(hdu, 1, 1.0, 0.0, "8A", NONULL, "ANNAME", "", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 2, 1.0, 0.0, "3D", NONULL, "STABXYZ", "METERS", NULL, NULL))
    return hduerr(hdu);
  sprintf(tform, "%dD", an->numorb);
  if(setbfield(hdu, 3, 1.0, 0.0, tform, NONULL, "ORBPARM", "", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 4, 1.0, 0.0, "1J", NONULL, "NOSTA", "", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 5, 1.0, 0.0, "1J", NONULL, "MNTSTA", "", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 6, 1.0, 0.0, "1E", NONULL, "STAXOF", "METERS", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 7, 1.0, 0.0, "1A", NONULL, "POLTYA", "", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 8, 1.0, 0.0, "1E", NONULL, "POLAA", "DEGREES", NULL, NULL))
    return hduerr(hdu);
  sprintf(tform, "%dE", an->nopcal);
  if(setbfield(hdu, 9, 1.0, 0.0, tform, NONULL, "POLCALA", "", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 10, 1.0, 0.0, "1A", NONULL, "POLTYB", "", NULL, NULL))
    return hduerr(hdu);
  if(setbfield(hdu, 11, 1.0, 0.0, "1E", NONULL, "POLAB", "DEGREES", NULL, NULL))
    return hduerr(hdu);
  sprintf(tform, "%dE", an->nopcal);
  if(setbfield(hdu, 12, 1.0, 0.0, tform, NONULL, "POLCALB", "", NULL, NULL))
    return hduerr(hdu);
/*
 * Add the initialized HDU to the FITS file.
 */
  if(add_Hdu(fits, hdu))
    return hduerr(hdu);
/*
 * Write the keywords describing the current sub-array. These were originally
 * copied from the input FITS file. If the source AN table was an ASCII
 * table then many of these parameters will not have been available
 * and this 
 *
 * Record the position of the center of the array.
 */
  if(wfltkey(fits, hdu, "ARRAYX", 0, an->arrayx,
	     "Array center X coord wrt Earth center (meters)") ||
     wfltkey(fits, hdu, "ARRAYY", 0, an->arrayy, "Array center Y coord.") ||
     wfltkey(fits, hdu, "ARRAYZ", 0, an->arrayz, "Array center Z coord."))
    return 1;
/*
 * Record the GST at time=0 on the ref. date.
 */
  if(wfltkey(fits, hdu, "GSTIA0", 0, an->gstia0,
	     "GST at time=0 on the ref. date (degrees)."))
    return 1;
/*
 * Record the Earth rotation rate.
 */
  if(wfltkey(fits, hdu, "DEGPDY", 0, an->degpdy,
	     "Earth rotation rate (deg/day)."))
    return 1;
/*
 * Record the reference frequency.
 */
  if(wfltkey(fits, hdu, "FREQ", 0, an->freq, "Ref. freq. of sub-array."))
    return 1;
/*
 * Record the reference date.
 */
  if(an->rdate[0] &&
     wstrkey(fits, hdu, "RDATE", 0, an->rdate, "Ref. date (DD/MM/YY)"))
    return 1;
/*
 * Polar positions X and Y.
 */
  if(wfltkey(fits, hdu, "POLARX", 0, an->polarx, "Polar X position.") ||
     wfltkey(fits, hdu, "POLARY", 0, an->polary, "Polar Y position."))
    return 1;
/*
 * Other date parameters.
 */
  if(wfltkey(fits, hdu, "UT1UTC", 0, an->ut1utc, "UT1-UTC (sec)") ||
     wfltkey(fits, hdu, "DATUTC", 0, an->datutc, "Data time-UTC (sec)"))
    return 1;
  if(an->timsys[0] &&
     wstrkey(fits, hdu, "TIMSYS", 0, an->timsys, "Time system"))
    return 1;
/*
 * Array name.
 */
  if(an->arrnam[0] &&
     wstrkey(fits, hdu, "ARRNAM", 0, an->arrnam, "Array name."))
    return 1;
/*
 * Number of orbital parameters, hard-wired to 0;
 */
  if(wintkey(fits, hdu, "NUMORB", 0, an->numorb, "Number of orbital parameters."))
    return 1;
/*
 * Number of polarization calibration constants.
 */
  if(wintkey(fits, hdu, "NOPCAL", 0, an->nopcal, "Number of polarization parameters."))
    return 1;
/*
 * Feed polarization parameterization.
 */
  if(an->poltype[0] && wstrkey(fits, hdu, "POLTYPE", 0, an->poltype,
				  "Feed polarization parameterization."))
    return 1;
/*
 * Write polarization P_REFANT and P_DIFFnn keywords if necessary.
 */
  if(wrt_p_refant(fits, hdu, sub))
    return 1;
/*
 * End the header.
 */
  if(end_header(fits, hdu))
    return 1;
/*
 * Write the data section of the table.
 */
  if(binandata(fits, (Thdu *)hdu, sub) || end_data(fits, hdu))
    return 1;
  return 0;
}

/*.......................................................................
 * Write the data portion of a binary AIPS AN antenna table.
 *
 * Input:
 *  fits       Fits *  The descriptor of the FITS file.
 *  thdu       Thdu *  The descriptor of the table HDU.
 *  sub    Subarray *  The sub-array to be recorded.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int binandata(Fits *fits, Thdu *thdu, Subarray *sub)
{
  Binan *an;    /* Binary AN table description */
  Bintel *tel;  /* Descriptor of the antenna being recorded */
  int irow;     /* The row in the table (telescope index + 1) */
/*
 * Get the container of input binary AN-table info.
 */
  an = sub->binan;
/*
 * Write one table row per antenna.
 */
  tel = an->bt;
  for(irow=1; irow<=sub->nstat; irow++,tel++) {
/*
 * Telescope name.
 */
    if(tel->anname[0] &&
       wcolumn(fits, thdu, 1, irow, DAT_CHR, 0, NULL, 0, 8, tel->anname)==0)
      return 1;
/*
 * Station X,Y,Z.
 */
    if(wcolumn(fits, thdu, 2, irow, DAT_DBL, 0, NULL, 0, 3, tel->stabxyz)==0)
      return 1;
/*
 * Orbital parameters.
 */
    if(an->numorb > 0 &&
       wcolumn(fits, thdu, 3, irow, DAT_DBL, 0, NULL, 0, an->numorb,
	       tel->orbparm)==0)
      return 1;
/*
 * Station number.
 */
    if(wcolumn(fits, thdu, 4, irow, DAT_INT, 0, NULL, 0, 1, &tel->nosta)==0)
      return 1;
/*
 * Mount type.
 */
    if(wcolumn(fits, thdu, 5, irow, DAT_INT, 0, NULL, 0, 1, &tel->mntsta)==0)
      return 1;
/*
 * Axis offset.
 */
    if(wcolumn(fits, thdu, 6, irow, DAT_DBL, 0, NULL, 0, 1, &tel->staxof)==0)
      return 1;
/*
 * Feed A feed poln.
 */
    if(wcolumn(fits, thdu, 7, irow, DAT_CHR, 0, NULL, 0, 1, &tel->poltya)==0)
      return 1;
/*
 * Feed A position angle.
 */
    if(wcolumn(fits, thdu, 8, irow, DAT_DBL, 0, NULL, 0, 1, &tel->polaa)==0)
      return 1;
/*
 * Feed A polarization cal parameters.
 */
    if(an->nopcal>0 && wcolumn(fits, thdu, 9, irow, DAT_DBL, 0, NULL, 0,
			       an->nopcal, tel->polcala)==0)
      return 1;
/*
 * Feed B feed poln.
 */
    if(wcolumn(fits, thdu, 10, irow, DAT_CHR, 0, NULL, 0, 1, &tel->poltyb)==0)
      return 1;
/*
 * Feed B position angle.
 */
    if(wcolumn(fits, thdu, 11, irow, DAT_DBL, 0, NULL, 0, 1, &tel->polab)==0)
      return 1;
/*
 * Feed B polarization cal parameters.
 */
    if(an->nopcal>0 && wcolumn(fits, thdu, 12, irow, DAT_DBL, 0, NULL, 0,
			       an->nopcal, tel->polcalb)==0)
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Return an 8 character '\0' terminated random parameter name for the
 * given UU,VV or WW parameter name. This includes appending a projection
 * code where relevant.
 *
 * Input:
 *  ob    Observation *   The descriptor of the observation.
 *  uvw          char *   UU,VV or WW.
 * Output:
 *  return       char *   The output random parameter name.
 */
static char *uvwname(Observation *ob, char *uvw)
{
  static char name[9];  /* The output string */
/*
 * Get the projection name.
 */
  char *proj = ob->proj==PRJ_SIN ? "" : Proj_name(ob->proj);
  sprintf(name, "%-4.4s%4.4s", uvw, proj);
/*
 * If a projection name was appended, convert interword spaces to hyphens.
 */
  if(*proj) {
    char *cptr;
    for(cptr=name; *cptr; cptr++) {
      if(*cptr == ' ')
	*cptr = '-';
    };
  };
  return name;
}

/*.......................................................................
 * If sub->p_refant is >= 0, write polarization P_REFANT and P_DIFFnn
 * keywords to the header of an antenna table.
 *
 * Input:
 *  fits       Fits *  The descriptor of the FITS file.
 *  hdu         Hdu *  The Hdu descriptor of the antenna table.
 *  sub    Subarray *  The sub-array associated with the antenna table.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int wrt_p_refant(Fits *fits, Hdu *hdu, Subarray *sub)
{
  int cif;
  if(sub->p_refant >= 0) {
/*
 * Record the reference antenna number.
 */
    if(wintkey(fits, hdu, "P_REFANT", 0, sub->p_refant, "Reference antenna"))
      return 1;
/*
 * Record the 'nif' R-L phase differences.
 */
    for(cif=0; cif<sub->nif; cif++) {
      if(wfltkey(fits, hdu, "P_DIFF", cif+1, sub->p_diff[cif],
		 "P_REFANT R-L phase difference"))
	return 1;
    };
  };
  return 0;
}
