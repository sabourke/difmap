#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "obs.h"
#include "units.h"
#include "vlbconst.h"
#include "vlbinv.h"
#include "vlbmath.h"
#include "mapmem.h"
#include "logio.h"

/*
 * Define the container of the Gridding Convolution Function.
 */
#define NGCF 301
typedef struct {
  float convfn[NGCF]; /* The sampled Gridding Convolution Function array */
  float tgtocg;       /* Converts from UV grid index to convfn[] index. */
} UVgcf;

static const int nmask=2;   /* The number of pixels on either side of a */
                            /*  given U,V into which to interpolate */

static int uvgrid(Observation *ob, MapBeam *mb, UVgcf *gcf, float uvmin,
	      float uvmax, float gauval, float gaurad, int dorad, float errpow,
	      int dounif, int domap);

static int uvbin(Observation *ob, MapBeam *mb, float binwid,
		 float uvmin, float uvmax);

static Bincell *getuvbin(UVbin *uvb, float uu, float vv);

static UVgcf *uvgcf(MapBeam *mb);

static void uv_limits(MapBeam *mb, float *ulimit, float *vlimit);
static int ispow2(int n);

/*.......................................................................
 * Fourier invert the residuals between the established model (after
 * establishing any tentative model) and the observed visibilities, to
 * yield a residual map, and/or fourier invert the UV plane sampling
 * to yield a dirty beam.
 *
 * Input:
 *  ob    Observation *  The UV data set. The tentative model will be
 *                       established before inverstion.
 *  uvmin       float    The UV radius (wavelengths) below which to ignore
 *                       data.
 *  uvmax       float    The UV radius (wavelengths) beyond which to ignore
 *                       data. If the largest of uvmin and uvmax is <= 0.0f
 *                       then the range will be unrestricted.
 *  gauval      float    The value of the weighting gaussian at UV radius
 *                       gaurad, between 0 and 1. If <=0 or >=1, no gaussian
 *                       taper is applied.
 *  gaurad      float    The radius (wavelengths) in the UV plane at which
 *                       the gaussian weighting function has value 'gauval'.
 *                       If <=0.0, no gaussian taper is applied.
 *  dorad         int    If true apply radial weighting in addition to
 *                       to gaussian weighting etc..
 *  errpow      float    If < 0.0 then the amplitude errors, raised to the
 *                       power 'errpow', will be used to scale the weights.
 *  binwid      float    For uniform weighting this specifies the width
 *                       of the square bin-size in UV pixels. Set to <=0
 *                       if uniform weighting is not required.
 * Input/Output:
 *  mb    MapBeam *      An intitialised map and beam container - see
 *                       function new_MapBeam() for details.
 * Output:
 *  return   int         0 - invert ended succesfully.
 *                       Anything else signifies an error.
 */
int uvinvert(Observation *ob, MapBeam *mb, float uvmin, float uvmax,
	     float gauval, float gaurad, int dorad, float errpow,
	     float binwid)
{
  UVgcf *gcf;    /* The gridding convolution function container */
  int old_if;    /* State of current IF to be restored on exit */
/*
 * Neither beam nor map has been requested - oops.
 */
  if(!mb->dobeam && !mb->domap) {
    lprintf(stderr, "uvinvert(): Neither beam nor a map requested.\n");
    return 1;
  };
/*
 * Check whether the observation is in an approriate state.
 */
  if(!ob_ready(ob, OB_SELECT, "uvinvert"))
    return 1;
/*
 * Store the state of the current IF.
 */
  old_if = get_cif_state(ob);
/*
 * Establish the tentative model if the map is to be computed.
 */
  if(mb->domap && mergemod(ob, 1))
    return 1;
/*
 * Inform the user.
 */
  lprintf(stdout, "Inverting %s%s%s\n", (mb->domap) ? "map " : "",
	 (mb->dobeam && mb->domap) ? "and " : "",
	 (mb->dobeam) ? "beam " : "");
/*
 * Mark the map as being dirty.
 */
  mb->ncmp = 0;
/*
 * Bin visibilities from all IFs, in preparation for uniform weighting in
 * uvgrid(). This function also checks the UV range in each IF, against
 * the grid size and must always be called. It will not actually bin the
 * data if binwid==0. This represents the case for natural weighting.
 */
  if(uvbin(ob, mb, binwid, uvmin, uvmax))
    return 1;
/*
 * Prepare the gridding interpolation function to be used in uvgrid(), and its
 * transform to be used in uvtrans().
 */
  gcf = uvgcf(mb);
  if(gcf==NULL)
    return 1;
/*
 * Grid the UV data into half of a conjugate symmetric array then
 * transform to the dirty map/beam.
 */
  if(mb->domap) {
    uvgrid(ob, mb, gcf, uvmin, uvmax, gauval, gaurad, dorad, errpow,
	   binwid>0, 1);
    uvtrans(mb, 1);
    mapstats(ob, mb);   /* Record the min/max valued pixels */
    mb->domap = 0;
  };
  if(mb->dobeam) {
    uvgrid(ob, mb, gcf, uvmin, uvmax, gauval, gaurad, dorad, errpow,
	   binwid>0, 0);
    uvtrans(mb, 0);
    mb->dobeam = 0;
  };
/*
 * Reinstate the original IF.
 */
  if(set_cif_state(ob, old_if))
    return 1;
  return 0;
}

/*.......................................................................
 * Interpolate the UV data points onto a UV grid of 'ngrid/2+1' by
 * 'ngrid' complex numbers, using a guassian interpolation function.
 * The result is in the form required to get back to the image plane
 * using realfft(), which returns the real valued map in a grid of
 * ngrid*ngrid floats. If you want the map to be centered then call
 * uvtrans() which phase shifts the array, then calls realfft().
 *
 * Input:
 *  ob Observation *  The observation to be gridded.
 *  mb     MapBeam *  The map and beam grid container.
 *  gcf      UVgcf *  The gridding convolution function descriptor.
 *  uvmin    float    The UV radius (wavelengths) below which to ignore
 *                    data.
 *  uvmax    float    The UV radius (wavelengths) beyond which to ignore
 *                    data. If the largest of uvmin and uvmax is <= 0.0f
 *                    then the range will be unrestricted.
 *  gauval   float    The value of the weighting gaussian at UV radius
 *                    gaurad, between 0 and 1. If <=0 or >=1, no gaussian
 *                    taper is applied.
 *  gaurad   float    The radius (wavelengths) in the UV plane at which
 *                    the gaussian weighting function has value 'gauval'.
 *                    If <=0.0, no gaussian taper is applied.
 *  dorad      int    If true apply radial weighting in addition to
 *                    to gaussian weighting etc..
 *  errpow   float    If < 0.0 then the amplitude errors, raised to the
 *                    power 'errpow', will be used to scale the weights.
 *  dounif     int    If true then uniform weighting will be performed.
 *  domap      int      0 - grid the dirty beam.
 *                      1 - grid the dirty map.
 * Output:
 *  mb->map or beam   An array of (ngrid+2)*ngrid floats to be treated as
 *                    the 2-D array of (ngrid/2+1)*ngrid real,imaginary
 *                    pairs, of the gridded UV points in one half of a
 *                    conjugate symettric array. NB. The 0,0 U,V point will
 *                    be located in array element 0,0.
 *  mb->e_bpa,        An estimate of the equivalent elliptical clean beam
 *  mb->e_bmin,       is recorded. All values are recorded in radians.
 *  mb->e_bmaj
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int uvgrid(Observation *ob, MapBeam *mb, UVgcf *gcf, float uvmin,
	      float uvmax, float gauval, float gaurad, int dorad, float errpow,
	      int dounif, int domap)
{
  float *uvmap;   /* Pointer to the map or beam array to be gridded */
  float *convfn;  /* The gridding convolution function array */
  int docut;      /* Flag whether uvmin and uvmax should be applied */
  int dotaper;    /* If true then apply a gaussian weighting taper */
  float gfac=0.0f;/* The reciprocal of the variance of the gaussian taper */
  float tgtocg;   /* Conversion factor between target grid pixels and */
                  /*  interpolation grid pixels. */
  float *cntr_ptr;/* Pointer to centre of 2D representation of uvmap */
  float *rptr;    /* Pointer to a complex element in the UV grid */
  float rval;     /* Interpolated real value of visibility */
  float ival;     /* Interpolated imaginary value of visibility */
  float wsum;     /* The sum of weights applied during gridding */
  float fv,fuv;   /* Value of interpolation function at centre of a pixel */
  float ulimit;   /* The largest |U| distance that can be accomodated */
  float vlimit;   /* The largest |V| distance that can be accomodated */
  int iu, iv;     /* U,V pixel coordinate in convolution grid */
  int nugrid;     /* Number of complex elements along U direction */
  int nvgrid;     /* The number of complex elements along V direction */
  int vinc;       /* Increment in floats to move up/down V axis */
  float *normptr; /* Pointer to U=0 at a given value of V=v */
  float *conjptr; /* Pointer to U=0 at V=-v (wrt normptr) */
  int cif;        /* The index of the IF being processed */
  int i;
/*
 * Create and initialize a container for the beam and noise estimation
 * sums and running means.
 */
  struct {
    float wsum;  /* Sum of un-interpolated gridding weights */
    float muu;   /* Mean of U.U */
    float mvv;   /* Mean of V.V */
    float muv;   /* Mean of U.V */
    float nsum;  /* Sum of grid weight / visibility weight */
  } bm = {0.0f,0.0f,0.0f,0.0f,0.0f};
/*
 * Get a pointer to the map or beam, the size of the grid and the number
 * of complex elements along the U axis.
 */
  uvmap = domap ? mb->map : mb->beam;
  nvgrid = mb->ny;
  nugrid = mb->nx/2+1;
/*
 * Set up for gridding - zero uvmap.
 */
  for(i=0; i<2*nugrid*nvgrid; i++)
    uvmap[i] = 0.0f;
/*
 * Enforce positivity on uvmin and uvmax.
 */
  if(uvmin < 0.0f)
    uvmin = 0.0f;
  if(uvmax < 0.0f)
    uvmax = 0.0f;
/*
 * Arrange that uvmin <= uvmax.
 */
  if(uvmin > uvmax) {float ftmp = uvmin; uvmin = uvmax; uvmax = ftmp;};
/*
 * Should we apply a cut-off in U and V?
 */
  docut = uvmax > 0.0f;
/*
 * Get the maximum U and V coordinates that can be nyquist sampled
 * using the current map pixel size.
 */
  uv_limits(mb, &ulimit, &vlimit);
/*
 * Record whether a gaussian taper was specified.
 */
  dotaper = gaurad > 0.0 && gauval > 0.0 && gauval < 1.0;
/*
 * Work out the -ve reciprocal of the variance of the gaussian taper.
 */
  if(dotaper)
    gfac = log(gauval)/gaurad/gaurad;
/*
 * Get the convolution function array, and the conversion factor between
 * pixels in the target grid and pixels in the convolution grid.
 */
  tgtocg = gcf->tgtocg;
  convfn = gcf->convfn;
/*
 * Get a pointer to the real part of the pixel U=0,V=N/2
 */
  cntr_ptr = uvmap + nvgrid * nugrid;
/*
 * Loop through all UV points in observation.
 */
  wsum = 0.0f;
/*
 * Loop through all sampled IFs.
 */
  for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
    float uvscale;  /* UV coordinate scale factor */
    Subarray *sub;  /* The descriptor of the sub-array being processed */
    int isub;       /* The index of sub in ob->sub[] */
/*
 * Get the next IF.
 */
    if(getIF(ob, cif))
      return 1;
/*
 * Get the multiplicative factor required to scale UVW light-second
 * distances to wavelength numbers at the frequency of the new IF.
 */
    uvscale = ob->stream.uvscale;
/*
 * Loop through sub-arrays of the new IF.
 */
    sub = ob->sub;
    for(isub=0; isub<ob->nsub; isub++,sub++) {
/*
 * Loop over integrations of the new sub-array.
 */
      Integration *integ = sub->integ;
      int ut;
      for(ut=0; ut<sub->ntime; ut++,integ++) {
/*
 * Loop over visibilities in the new integration.
 */
	Visibility *vis = integ->vis;
	int base;
	for(base=0; base<sub->nbase; base++,vis++) {
          float uu = vis->u * uvscale;       /* U distance (wavelengths) */
          float vv = vis->v * uvscale;       /* V distance (wavelenghts) */
          float uvrad = sqrt(uu*uu + vv*vv); /* Radial distance in UV plane */
/*
 * Only grid usable visibilities.
 */
	  if(!vis->bad && !(docut && (uvrad < uvmin || uvrad > uvmax)) &&
	     fabs(uu) <= ulimit && fabs(vv) <= vlimit) {
	    float ufrc = uu / mb->uinc;   /* Decimal pixel position */
	    float vfrc = vv / mb->vinc;
	    int upix = fnint(ufrc);       /* Integer pixel position */
	    int vpix = fnint(vfrc);
	    float uvrval;                 /* Real part of visibility */
	    float uvival;                 /* Imaginary part of visibility */
	    float weight=1.0f;            /* Weight to apply to visibility */
/*
 * Work out the weight to assign to the new visibility.
 */
	    if(dotaper)
	      weight *= exp(gfac * uvrad*uvrad);  /* Gaussian taper. */
/*
 * Radial weighting.
 */
	    if(dorad)
	      weight *= uvrad;
/*
 * Amplitude uncertainty weighting - include special cases for common
 * values.
 */
	    if(errpow < -0.001) {
	      float power = -errpow/2.0f;
	      float wt = fabs(vis->wt);
	      if(power==1.0f)
		weight *= wt;             /* vis->wt is the correct value */
	      else if(power==0.5f)
		weight *= sqrt(wt);       /* sqrt() is faster than pow() */
	      else
		weight *= pow(wt, power); /* General case */
	    };
/*
 * Uniform weighting?
 */
	    if(dounif) {
	      Bincell *bc = getuvbin(mb->bin, uu, vv);
	      if(bc && *bc>0)
		weight /= *bc;
	    };
/*
 * Accumulate the weighted running means used to estimate the clean-beam.
 * Use of running means is essential since the numbers being added are
 * very large.
 */
	    if(!domap) {
	      float runwt = weight / (bm.wsum += weight);
	      bm.muu += runwt * (uu*uu - bm.muu);
	      bm.mvv += runwt * (vv*vv - bm.mvv);
	      bm.muv += runwt * (uu*vv - bm.muv);
/*
 * Accumulate weight sum used together with bm.wsum to calculate the
 * estimated noise.
 */
	      bm.nsum += weight * weight / vis->wt;
	    };
/*
 * Turn the data value into a complex form.
 */
	    if(domap) {
	      uvrval = vis->amp * cos(vis->phs) -
		vis->modamp * cos(vis->modphs);
	      uvival = vis->amp * sin(vis->phs) -
		vis->modamp * sin(vis->modphs);
	    } else {
	      uvrval = 1.0f;    /* Beam */
	      uvival = 0.0f;
	    };
/*
 * Convolve the 2*nmask+1 square array of points around upix and vpix
 * with the interpolation function. In the full conjugate symmetric array
 * each point is mirrored by its conjugate value on the opposite side of the
 * origin. In the half-array that we are building only +ve U values
 * are stored (since the symmetry makes it redundant to store the other side).
 * Where a point is located in the -ve U part of the plane, its conjugate
 * mirror image will be inserted instead.
 *
 * Loop through the interpolation area.
 */
	    for(iv = vpix-nmask; iv<=vpix+nmask; iv++) {
/*
 * Determine the value of the interpolation function along V at this pixel.
 */
	      fv = weight * convfn[(int) (tgtocg*fabs(iv-vfrc)+0.5f)];
/*
 * Determine the increment in floats to move from v=N/2 to v=vpix+iv.
 * The same increment with the opposite sign will take us to v=-N/2, (except
 * when v=0 [see below]) hence the choice of U=0,V=N/2 as the reference point.
 */
	      vinc = nugrid*(iv+iv+((iv<0)?nvgrid:-nvgrid));
/*
 * Determine pointers to U=0,V=iv and U=0,V=-iv.
 */
	      normptr = cntr_ptr + vinc;
	      conjptr = cntr_ptr + ((iv) ? -vinc:vinc);	      
	      for(iu = upix-nmask; iu<=upix+nmask; iu++) {
/*
 * Combine the interpolation functions along U and V.
 */
		wsum += (fuv = fv * convfn[(int) (tgtocg*fabs(iu-ufrc)+0.5f)]);
/*
 * Calculate the real and imaginary parts of the interpolated
 * and weighted UV data value.
 */
		rval = uvrval*fuv;
		ival = uvival*fuv;
/*
 * Pixel iu,iv may be inside the array or in the non-existent
 * conjugate other half of the array. If it is in the latter
 * then we should put it at its conjugate symmetric position in
 * the array - this also means that the gridded data value should be
 * conjugated.
 */
		if(iu <= 0) {
		  rptr = conjptr-iu-iu; /* Pointer to conjugate element */
		  *rptr += rval;
		  *(rptr+1) -= ival;
		};
		if(iu >= 0) {
		  rptr = normptr+iu+iu; /* Pointer to complex element */
		  *rptr += rval;
		  *(rptr+1) += ival;
		};
	      };
	    };
	  };
	};
      };
    };
  };
/*
 * If a zero spacing flux has been specified convolve it in separately here
 * using the same algorithm as above (without the imaginary parts).
 * Note that the zero baseline flux has zero weight if radial weighting
 * has been selected, and should then be ignored.
 */
  if(ob->uvzero.wt > 0.0f && !dorad) {
    float weight = 1.0f;
    float uvrval = domap ? (ob->uvzero.amp - ob->uvzero.modamp) : 1.0f;
/*
 * Apply amplitude uncertainty weighting?
 */
    if(errpow < -0.001)
      weight *= pow(ob->uvzero.wt, -errpow/2.0f); /* General case */
/*
 * Uniform weighting?
 */
    if(dounif) {
      Bincell *bc = getuvbin(mb->bin, 0.0f, 0.0f);
      if(bc && *bc>0)
	weight /= *bc;
    };
/*
 * Convolve over the nmask pixels either side of cntr_ptr.
 */
    for(iv = -nmask; iv<=nmask; iv++) {
      fv = weight * convfn[(int) (tgtocg*fabs(iv)+0.5f)];
/*
 * Address increment to U=0,V=iv.
 */
      vinc = nugrid*(iv+iv+((iv<0)?nvgrid:-nvgrid));
/*
 * Determine pointers to U=0,V=iv and U=0,V=-iv.
 */
      normptr = cntr_ptr + vinc;
      conjptr = cntr_ptr + ((iv) ? -vinc:vinc);
      for(iu = -nmask; iu<=nmask; iu++) {
	wsum += (fuv = fv * convfn[(int) (tgtocg*fabs(iu)+0.5f)]);
	rval = uvrval * fuv;
	if(iu <= 0)
	  *(conjptr-iu-iu) += rval;   /* Conjugate element */
	if(iu >= 0)
	  *(normptr+iu+iu) += rval;   /* Sampled element */
      };
    };
  };
/*
 * No data gridded?
 */
  if(wsum<=0.0f || (!domap && bm.wsum<=0.0f)) {
    lprintf(stderr, "uvgrid: No data in UV range.\n");
    return 1;
  };
/*
 * Finally - divide the UV grid by the sum of weights. Prescale
 * wsum by 2 to take into account the fact that every point appears twice
 * in the UV plane.
 */
  wsum *= 2.0f;
  rptr=uvmap;
  for(iv=0; iv<nvgrid; iv++) {
    for(iu=0;iu<2*nugrid;iu++)
      *(rptr++) /= wsum;
  };
/*
 * Work out the estimate of the size of the clean beam.
 * The technique used was developed by Tim Pearson, and I don't fully
 * understand it. It depends on the property of fourier transforms
 * that relates the 2nd moment in the UV plane to the curvature at the center 
 * of the beam in the image plane. An empirical fudge factor is used to
 * extrapolate the extents of the beam at HWHM.
 */
  if(!domap) {
    const float fudge=0.7f; /* Empirical fudge factor of TJP's algorithm */
    float ftmp = sqrt((bm.muu-bm.mvv)*(bm.muu-bm.mvv) + 4.0*bm.muv*bm.muv);
/*
 * First the position angle of the equivalent elliptical gaussian distribution.
 */
    mb->e_bpa = -0.5*atan2(2.0*bm.muv, bm.muu - bm.mvv);
/*
 * Then the equivalent elliptical beam widths in radians.
 */
    mb->e_bmin = fudge/(sqrt(2.0*(bm.muu+bm.mvv) + 2.0*ftmp));
    mb->e_bmaj = fudge/(sqrt(2.0*(bm.muu+bm.mvv) - 2.0*ftmp));
    lprintf(stdout,
         "Estimated beam: bmin=%.4g %s, bmaj=%.4g %s, bpa=%.4g degrees\n",
	    radtoxy(mb->e_bmin), mapunits(U_NAME),
	    radtoxy(mb->e_bmaj), mapunits(U_NAME),
	    mb->e_bpa * rtod);
/*
 * Determine the estimated map noise.
 */
    mb->noise = sqrt(bm.nsum / bm.wsum / bm.wsum);
/*
 * Display the estimated noise.
 */
    lprintf(stdout, "Estimated noise=%g mJy/beam.\n", mb->noise * 1.0e+3);
  };
  return 0;
}

/*.......................................................................
 * Accumulate visibility counts for all IFs, binned in U and V for use in
 * uniform weighting. This function also checks the UV range in each IF
 * against the chosen UV grid size, so it must be called before calling
 * uvgrid(), regardless of whether uniform weighting is desired.
 *
 * The bin array in mb->uvbin has size mb->nx/4 x mb->ny/2 and is used
 * to bin the +ve U half of the conjugate-symmetric UV plane over the
 * UV range: U = 0 to |UMAX|/2, V = -|VMAX|/2 to |VMAX|/2. The extra factor
 * of a half is due to nyquist constraints that require that the UV plane
 * be no more than half sampled, to avoid undersampling the image plane.
 * If visibilities within the passed uvmin and uvmax values are discovered
 * this function will abort after informing the user.
 *
 * The actual number of bin elements used is further cut if binwid is
 * > 1, to accomodate larger UV bin widths.
 * 
 * Input:
 *  ob  Observation *  The observation to be uniform weighted.
 *  mb      MapBeam *  The descriptor containing the uniform weighting
 *                     array and other gridding parameters.
 *  binwid    float    The size of the uniform bin in grid pixels numbers.
 *                     If binwid==0, 1 will be placed in the U=0,V=0 bin
 *                     and uvb->utopix and uvb->vtopix will be set to 0.0
 *                     this represents natural weighting.
 *  uvmin     float    The minimum UV radius to take visibilities from.
 *  uvmax     float    The maximum UV radius to take visibilities from.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int uvbin(Observation *ob, MapBeam *mb, float binwid,
		 float uvmin, float uvmax)
{
  UVbin *uvb;    /* The UV bin descriptor, mb->bin. */
  Bincell *bc;   /* Pointer to a cell in the bin array uvb->bin[] */
  float ulimit;  /* The largest |U| distance that can be accomodated */
  float vlimit;  /* The largest |V| distance that can be accomodated */
  float umax;    /* The maximum u distance > ulimit */
  float vmax;    /* The maximum v distance > vlimit */
  int dounif;    /* True if uniform weighting is required */
  int docut;     /* Flag whether uvmin and uvmax should be applied */
  int cif;       /* The index of the IF being processed */
  int isub;      /* The index of the sub-array being processed */
  int ut;        /* The index of the integration being processed. */
  int base;      /* The index of the baseline being processed. */
  int ngood;     /* The number of good visibilities */
  int nused;     /* The number of ngood that were within the specified uv */
                 /*  radius range and within specified u and v maxima */
  int nbadr;     /* The number of ngood that were rejected due to being */
                 /*  outside of the specified uv radius range. */
  int nbaduv;    /* The number of ngood that were within the specified */
                 /*  uv radius range, but beyond the specified U or V maxima. */
  int i;
/*
 * Get the UV bin descriptor created by new_Mapmem().
 */
  uvb = mb->bin;
/*
 * Check the requested bin size.
 */
  if(binwid<0 || binwid >= uvb->nu || binwid >= uvb->nv) {
    lprintf(stderr, "Uniform bin width (%g) out of permissible range.\n",
	    binwid);
    return 1;
  };
/*
 * Is uniform weighting required?
 */
  dounif = binwid > 0.0;
/*
 * Uniform weigthing is constrained by the size of the weights array
 * to binwid >= 1.0.
 */
  if(dounif && binwid < 1.0) {
    lprintf(stderr, "Uniform bin width adjusted to minimum of 1.0.\n");
    binwid = 1.0;
  };
/*
 * Enforce positivity on uvmin and uvmax.
 */
  if(uvmin < 0.0f)
    uvmin = 0.0f;
  if(uvmax < 0.0f)
    uvmax = 0.0f;
/*
 * Arrange that uvmin <= uvmax.
 */
  if(uvmin > uvmax) {float ftmp = uvmin; uvmin = uvmax; uvmax = ftmp;};
/*
 * Should we apply a cut-off in UV radius?
 */
  docut = uvmax > 0.0f;
/*
 * Determine the conversion factor between U and V (wavelength) coords
 * to bin array indexes. This is used in subsequent calls to getuvbin().
 */
  uvb->utopix = dounif ? 1.0f / mb->uinc / binwid : 0.0f;
  uvb->vtopix = dounif ? 1.0f / mb->vinc / binwid : 0.0f;
/*
 * Get the maximum U and V coordinates that can be nyquist sampled
 * using the current map pixel size.
 */
  uv_limits(mb, &ulimit, &vlimit);
/*
 * Initialize counters.
 */
  ngood = 0;
  nused = 0;
  nbadr = 0;
  nbaduv = 0;
/*
 * Initialize maxima.
 */
  umax = vmax = 0.0f;
/*
 * Zero the work array.
 */
  bc = uvb->bins;
  for(i=0; i<uvb->nbin; i++)
    *bc++ = 0;
/*
 * Loop through all sampled IFs.
 */
  for(cif=0; (cif=nextIF(ob, cif, 1, 1)) >= 0; cif++) {
    float uvscale;   /* Scale factor from UV coords to wavelengths */
    Subarray *sub;   /* The current sub-array */
/*
 * Get the next IF to be visited.
 */
    if(getIF(ob, cif))
      return 1;
/*
 * Get the conversion factor between UV coords and wavelengths in the
 * new IF.
 */
    uvscale = ob->stream.uvscale;
/*
 * For the U,V coordinate of each visibility locate its equivalent pixel
 * in the UV array and add 1.0 to that location.
 */
    sub = ob->sub;
    for(isub=0; isub<ob->nsub; isub++, sub++) {
      Integration *integ = sub->integ;
      for(ut=0; ut<sub->ntime; ut++,integ++) {
	Visibility *vis = integ->vis;
	for(base=0; base<sub->nbase; base++,vis++) {
/*
 * Ignore flagged visibilities.
 */
	  if(!vis->bad) {
	    float uu = vis->u * uvscale;
	    float vv = vis->v * uvscale;
	    float fabs_uu = fabs(uu);
	    float fabs_vv = fabs(vv);
	    float uvrad = sqrt(uu*uu + vv*vv);
/*
 * Count unflagged visibilities.
 */
	    ngood++;
/*
 * Reject and count any visibilities that fall outside the specified
 * uv radius range.
 */
	    if(docut && (uvrad < uvmin || uvrad > uvmax)) {
	      nbadr++;
/*
 * Reject visibilities that exceed the allowed range of U and V.
 */
	    } else if(fabs_uu > ulimit || fabs_vv > vlimit) {
	      nbaduv++;
/*
 * Record the worst overflow of the allowed UV range.
 */
	      if(fabs_uu > umax)
		umax = fabs_uu;
	      if(fabs_vv > vmax)
		vmax = fabs_vv;
/*
 * Process acceptable visibilities.
 */
	    } else {
	      nused++;
/*
 * Collect uniform weighting bin counts?
 */
	      if(dounif) {
		bc = getuvbin(uvb, uu, vv);
		if(bc) (*bc)++;
/*
 * If the visibility is in the U=0 bin then its conjugate mirrored
 * point will also be in the bin array, on the other side of V=0.
 * Count it as well.
 */
		if(fnint(fabs_uu * uvb->utopix)==0) {
		  bc = getuvbin(uvb, uu, -vv);
		  if(bc) (*bc)++;
		};
	      };
	    };
	  };
	};
      };
    };
  };
/*
 * Add a uniform bin entry to account for the optional zero-spacing flux
 * or for natural weighting.
 */
  bc = getuvbin(uvb, 0.0f, 0.0f);
  if(bc) (*bc)++;
/*
 * Are all of the visibilities flagged?
 */
  if(ngood <= 0) {
    lprintf(stderr, "There are no unflagged visibilities to be inverted.\n");
    return 1;
  };
/*
 * Report the number of visibilities that were excluded because of the
 * chosen uv radius range.
 */
  if(nbadr) {
    lprintf(stderr,
	    "Your chosen uvrange limits excluded %.2g%% of the data.\n",
	    100.0 * nbadr / (double)ngood);
  };
/*
 * Report the number of visibilities that were excluded because the
 * current map cell size undersamples those visibilities.
 */
  if(nbaduv) {
    lprintf(stderr,
	    "Your choice of large map pixels excluded %s%.3g%% of the data.\n",
	    nbadr ? "a further ":"", 100.0 * nbaduv / (double)ngood);
    if(umax > ulimit) {
      lprintf(stderr,
	      " The x-axis pixel size should ideally be below %.4g %s\n",
	      radtoxy(ulimit/umax * mb->xinc), mapunits(U_TLAB));
    };
    if(vmax > vlimit) {
      lprintf(stderr,
	      " The y-axis pixel size should ideally be below %.4g %s\n",
	      radtoxy(vlimit/vmax * mb->xinc), mapunits(U_TLAB));
    };
  };
/*
 * Did non of the unflagged visibilities lie within the current uv ranges?
 */
  if(nused <= 0) {
    lprintf(stderr, "No visibilities were available for creating a map.\n");
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Return a pointer to the UV bin corresponding to a given U and V position.
 * uvbin() must have been called before this function, to initialize
 * binning parameters.
 *
 * Input:
 *  uvb     UVbin *  The UV bin array descriptor.
 *                   This must NOT be changed between calling uvbin() and
 *                   uvgrid().
 *  uu      float    The U coordinate to acquire a bin for (wavelengths).
 *  vv      float    The V coordinate to acquire a bin for (wavelengths).
 * Output:
 *  return Bincell * The pointer to the bin cell, or NULL on error.
 */
static Bincell *getuvbin(UVbin *uvb, float uu, float vv)
{
  long binpix;  /* The 1-D array index into the 2-D bin array */
/*
 * If the point is in the unmapped -ve U part of the UV plane,
 * switch to the mapped conjugate-symmetric pixel.
 */
  if(uu >= 0) {
    uu = -uu;
    vv = -vv;
  };
/*
 * Determine the position in the bin array wrt its U=0,V=0 origin.
 */
  binpix = uvb->nu * (uvb->nv/2 + floor(vv * uvb->vtopix + 0.5)) +
    floor(uu * uvb->utopix + 0.5);
/*
 * Return the corresponding address.
 */
  return binpix >= 0 && binpix < uvb->nbin ? uvb->bins + binpix : NULL;
}

/*.......................................................................
 * Calculate the interpolation function required to convolve visibilities
 * onto the UV grid - usually called the Gridding Convolution Function (gcf).
 * This will be returned for use in uvgrid(). Also record its normalized
 * Fourier transform in mb->rxft[] (X-axis), and mb->ryft[] (Y-axis) for
 * use in uvtrans().
 *
 * Input:
 *  mb    MapBeam *  The container of information pertaining to the
 *                   UV grid. Parameters used are nx and ny. Also
 *                   mb->rxft[] and mb->ryft[] are initialized.
 * Output:
 *  return  UVgcf *  A pointer to the static internal container of the
 *                   convolution array, for use in uvgrid().
 *                   Or NULL on error.
 */
static UVgcf *uvgcf(MapBeam *mb)
{
  static UVgcf gcf;       /* Convolution function container */
  const float hwhm=0.7f;  /* The HWHM of the convolution gaussian in */
                          /* multiples of elements on the target grid. */
  float recvar;   /* The reciprocal of twice the gaussian function variance */
  float cghwhm;   /* HWHM of gaussian in gcf array pixels */
  float peak;     /* The peak in the un-normalized transorm of the gcf */
  float *rptr;    /* Pointer into mb->rxft[] or mb->ryft[] */
  int i;
/*
 * Check arguments.
 */
  if(mb==NULL) {
    lprintf(stderr, "uvgcf: NULL MapBeam descriptor intercepted.\n");
    return NULL;
  };
/*
 * Get a conversion factor between pixels in the target grid and
 * pixels in the convolution grid. The convolution is performed
 * to the 'nmask' pixels on either side of the closest pixel to the
 * UV point. The convolution grid thus corresponds to nmask+0.5
 * pixels on one side of the centre of the closest pixel. The
 * -1 in (NGCF-1) is a precaution, in that insufficient float precision
 * in indexing the gcf[] array while convolving could just
 * possibly result in accessing one element further than intended.
 */
  gcf.tgtocg = (NGCF-1)/(nmask+0.5);
/*
 * Convert the half-wid-half-max from a multiple of UV grid pixels
 * to a multiple of convolution grid pixels.
 */
  cghwhm = gcf.tgtocg * hwhm;
/*
 * Convert the half-width-half-maximum to the reciprocal of twice the
 * equivalent gaussian variance.
 */
  recvar = log(2.0)/cghwhm/cghwhm;
/*
 * Calculate the gaussian convolution function.
 */
  for(i=0; i<NGCF; i++)
    gcf.convfn[i] = exp(-recvar*i*i);
/*
 * Cosine transform the convolution function for both the X and Y
 * map/beam axes.
 */
  costran(gcf.convfn, NGCF-1, nmask+0.5, mb->rxft, mb->nx);
  costran(gcf.convfn, NGCF-1, nmask+0.5, mb->ryft, mb->ny);
/*
 * Take the normalized reciprocals of the FT of the convolution function.
 * The result can then be used in uvtrans() to deconvolve the convolution
 * function.
 */
  peak = mb->rxft[mb->nx/2];  /* Central peak value */
  rptr = mb->rxft;
  for(i=0; i<mb->nx; i++,rptr++)
    *rptr = peak / *rptr;
/*
 * Now the Y-axis transform.
 */
  peak = mb->ryft[mb->ny/2];  /* Central peak value */
  rptr = mb->ryft;
  for(i=0; i<mb->ny; i++,rptr++)
    *rptr = peak / *rptr;
/*
 * Return the container of the convolution function grid.
 */
  return &gcf;
}

/*.......................................................................
 * Return the maximum U and V coordinates that can be nyquist sampled
 * using the current map pixel size.
 *
 * Nyquist sampling of the image plane requires that only half of the
 * grid on either side of U=0 and V=0 be sampled.
 *
 * Input:
 *  mb     MapBeam *   The map and beam container object.
 * Input/Output:
 *  ulimit   float *   The maximum U coordinate that can be sampled
 *                     will be assigned to *ulimit.
 *  vlimit   float *   The maximum V coordinate that can be sampled
 *                     will be assigned to *vlimit.
 */
static void uv_limits(MapBeam *mb, float *ulimit, float *vlimit)
{
  *ulimit = mb->uinc * (mb->nx / 4 - nmask);
  *vlimit = mb->vinc * (mb->ny / 4 - nmask);
}

/*.......................................................................
 * Return the maximum pixel sizes at which all visibilities can be used
 * during gridding.
 *
 * Input:
 *  ob     Observation *  The observation to be characterized.
 *  uvmin        float    The minimum UV radius to take visibilities from.
 *  uvmax        float    The maximum UV radius to take visibilities from.
 *  nx, ny         int    The desired number of pixels in the map (including
 *                        the unseen margins).
 * Input/Output:
 *  xmax, ymax   float *  The maximum pixel sizes usable without loss of
 *                        data (radians).
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
int optimal_pixel_size(Observation *ob, float uvmin, float uvmax,
		       int nx, int ny, float *xmax, float *ymax)
{
  UVrange *uvr;   /* The range of data available */
/*
 * Check whether the observation is in an approriate state.
 */
  if(!ob_ready(ob, OB_SELECT, "optimal_pixel_size"))
    return 1;
/*
 * Check the arguments.
 */
  if(nx < 1 || !ispow2(nx) || ny < 1 || !ispow2(ny)) {
    lprintf(stderr, "Invalid non-power of 2 number of pixels.\n");
    return 1;
  };
/*
 * Get the range of data available.
 */
  uvr = uvrange(ob, 1, 0, uvmin, uvmax);
  if(!uvr)
    return 1;
/*
 * Make sure that there is some data.
 */
  if(uvr->umax <= 0.0 || uvr->vmax <= 0.0) {
    lprintf(stderr, "No data are in the current UV range.\n");
    return 1;
  };
/*
 * Work out the maximum pixel sizes along the x and y axes.
 */
  if(xmax)
    *xmax = (nx / 4 - nmask) / (uvr->umax * nx);
  if(ymax)
    *ymax = (ny / 4 - nmask) / (uvr->vmax * nx);
  return 0;
}

/*.......................................................................
 * Determine whether a given integer is a power of 2. This method is a
 * lot faster and more robust than using logs etc..
 *
 * Input:
 *  n        int   The grid size.
 * Output:
 *  return   int   0 - n is not a power of 2.
 *                 1 - n is a power of 2.
 */
static int ispow2(int n)
{
/*
 * Negative numbers are not powers of 2.
 */
  if(n<=0)
    return 0;
/*
 * Divide by 2 until the remainder of dividing by two is non-zero, or
 * 1 is reached. 1 is only reached if the number is a power of two.
 */
  while(n>1 && n%2==0)
    n >>= 1;
  return n==1;
}

