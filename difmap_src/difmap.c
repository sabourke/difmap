/*
 * This file contains user functions concerned with VLBI difference mapping
 * and user accessible variables used by these functions.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "logio.h"
#include "sphere.h"
#include "scrfil.h"
#include "helpdir.h"

#include "obs.h"
#include "vlbinv.h"
#include "mapmem.h"
#include "vlbconst.h"
#include "vlbutil.h"
#include "vlbmath.h"
#include "mapwin.h"
#include "mapcln.h"
#include "mapres.h"
#include "maplot.h"
#include "winmod.h"
#include "obwin.h"
#include "slfcal.h"
#include "wmap.h"
#include "telspec.h"
#include "visplot.h"
#include "scans.h"
#include "units.h"
#include "enumpar.h"
#include "baselist.h"
#include "pollist.h"
#include "specplot.h"
#include "cpgplot.h"
#include "modeltab.h"
#include "markerlist.h"
#include "visstat.h"
#include "planet.h"
#include "pb.h"
#include "mapcor.h"

extern char *date_str(void);

/*
 * Declare variables that are to be aliased as user variables below. Only
 * float, integer, char and logical variables are supported.
 * NB. character strings must NOT be initialised here unless they are marked
 * as R_ONLY parameters. This is to allow variable length strings where
 * the previous string is often free'd first on the assumption that the
 * memory for the string was allocated using malloc(), not by the
 * compiler).
 */

/* Struct to contain parameters for invert */

static struct {
  float uvhwhm;   /* HWHM of UV interpolation function (pixels) */
  float uvmin;    /* UV min radius cutoff (wavelengths) */
  float uvmax;    /* UV min radius cutoff (wavelengths) */
  float gauval;   /* Value of gaussian taper at 'gaurad' */
  float gaurad;   /* UV radius of Gaussian taper (wavelengths) */
  float errpow;   /* Exponent applied to errors for weighting */
  float uvbin;    /* UV bin width for uniform weighting */
  char dorad;     /* If true, do radial weighting */
} invpar, invdef={0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,2.0f,0};

/* Struct to contain self-cal parameters */

static struct {
  float gauval;   /* Value of gaussian taper at 'gaurad' */
  float gaurad;   /* UV radius of Gaussian taper (wavelengths) */
  float maxamp;   /* Max allowable amplitude correction */
  float maxphs;   /* Max allowable phase correction */
  int p_mintel;   /* Minimum number of telescopes for phase solution */
  int a_mintel;   /* Minimum number of telescopes for amplitude solution */
  int doflag;     /* If true, flag un-correctable baslines */
} slfpar, slfdef={0.0f,0.0f,0.0f,0.0f,3,4,1};

/* Struct to contain parameters for maplot */

static struct {
  Ctable *ctab;      /* Color table. */
  char docont;       /* If true then plot contours */
  float cmul;        /* Contour multiplier */
  float box[4];      /* Displayed area (radians) */
  MaplotBeam mpb;    /* Mapplot clean beam parameters */
  MaplotVect vect;   /* Polarization vector display parameters */
} mappar = {
  NULL, 1, 0.0f, {0.0f,0.0f,0.0f,0.0f}, {0.0f,0.0f,0.01f,0.3f},
  {0.0, 0.0, 0.0, 1, 1},
};

static int make_polmap(int docln);
static int polmap_error(Stokes pol);

/* Struct to contain parameters for mapres */

static struct {
  float bmin;    /* Beam semi-minor axis (radians) */
  float bmaj;    /* Beam semi-major axis (radians) */
  float bpa;     /* Major axis position angle (degrees) */
  float e_bmin;  /* Estimate of bmin from the last 'invert' (radians) */
  float e_bmaj;  /* Estimate of bmaj from the last 'invert' (radian) */
  float e_bpa;   /* Estimate of bpa from the last 'invert' (degrees) */
  int doauto;    /* Use estimate of beam size from last invert if true */
} respar, resdef={0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1};

/* Struct containing parameters for clean */

static struct {
  int niter;     /* Default max number of iterations */
  float gain;    /* Default CLEAN gain */
  float cutoff;  /* Default residual cut off */
} clnpar = {
  100, 0.05f, 0.0f
};

/* The following variable is set by the multi_model command */

static int multi_model_mode = 0;

static Descriptor mb_levs = {'f' , '1' ,RWD    ,1, {1,1,1}, NULL};
static Descriptor mb_map  = {'f' , '2' ,NO_DEL ,1, {1,1,1}, NULL};
static Descriptor mb_beam = {'f' , '2' ,NO_DEL ,1, {1,1,1}, NULL};
static Descriptor vflags  = {'c' , '0' ,NO_DEL ,1, {1,1,1}, NULL};
static Descriptor rflags  = {'c' , '0' ,NO_DEL ,1, {1,1,1}, NULL};
static Descriptor pflags  = {'c' , '0' ,NO_DEL ,1, {1,1,1}, NULL};
static Descriptor tflags  = {'c' , '0' ,NO_DEL ,1, {1,1,1}, NULL};
static Descriptor uflags  = {'c' , '0' ,NO_DEL ,1, {1,1,1}, NULL};

static Descriptor dmapv_type[] = {
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.uvhwhm },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.uvmin },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.uvmax },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.gauval },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.gaurad },
   {'l' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.dorad },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.errpow },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &invpar.uvbin },
   {'f' , '0' ,RWD    ,1, {1,1,1}, NULL },
   {'f' , '0' ,RWD    ,1, {1,1,1}, NULL },
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &mappar.cmul },
   {'l' , '0' ,RWD    ,1, {1,1,1}, NULL },
   {'l' , '0' ,NO_DEL ,1, {1,1,1}, &mappar.docont },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &respar.bmin },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &respar.bmaj },
   {'f' , '0' ,R_ONLY ,1, {1,1,1}, &respar.bpa },
   {'i' , '0' ,NO_DEL ,1, {1,1,1}, &clnpar.niter },
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &clnpar.gain },
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &clnpar.cutoff },
   {'D' , '1' ,NO_DEL ,1, {1,1,1}, &mb_levs },
   {'D' , '1' ,NO_DEL ,1, {1,1,1}, &mb_map },
   {'D' , '1' ,NO_DEL ,1, {1,1,1}, &mb_beam },
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &vflags },
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &rflags },
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &pflags },
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &tflags },
   {'D' , '0' ,NO_DEL ,1, {1,1,1}, &uflags },
};

/*
  In the same order as the above array of types, define the array
  of user names for the arrays.
*/

static char *dmapv_name[] = {
  "hwhm",
  "uvmin",
  "uvmax",
  "gauval",
  "gaurad",
  "dorad",
  "errpow",
  "uvbin",
  "gmin",
  "gmax",
  "cmul",
  "dogrey",
  "docont",
  "bmin",
  "bmaj",
  "bpa",
  "niter",
  "gain",
  "cutoff",
  "levs",
  "map",
  "beam",
  "vflags",
  "rflags",
  "pflags",
  "tflags",
  "uflags",
};

/*
 * Declare the user functions here.
 */

static Template(newob_fn);
static Template(mapsize_fn);
static Template(invert_fn);
static Template(uvtaper_fn);
static Template(uvrange_fn);
static Template(uvwgt_fn);
static Template(uvzero_fn);
static Template(maplot_fn);
static Template(clean_fn);
static Template(restore_fn);
static Template(wbeam_fn);
static Template(wmap_fn);
static Template(wdmap_fn);
static Template(wobs_fn);
static Template(wmodel_fn);
static Template(rmodel_fn);
static Template(gscal_fn);
static Template(keep_fn);
static Template(clrmod_fn);
static Template(shift_fn);
static Template(unshift_fn);
static Template(uvav_fn);
static Template(head_fn);
static Template(uncal_fn);
static Template(corpl_fn);
static Template(tname_fn);
static Template(ntel_fn);
static Template(bname_fn);
static Template(nbase_fn);
static Template(nsub_fn);
static Template(nif_fn);
static Template(nchan_fn);
static Template(addwin_fn);
static Template(delwin_fn);
static Template(winmod_fn);
static Template(startmod_fn);
static Template(uvrad_fn);
static Template(uvprj_fn);
static Template(self_fn);
static Template(staper_fn);
static Template(slims_fn);
static Template(rwins_fn);
static Template(wwins_fn);
static Template(vplot_fn);
static Template(uvplt_fn);
static Template(timpl_fn);
static Template(resof_fn);
static Template(unoff_fn);
static Template(save_fn);
static Template(get_fn);
static Template(loglev_fn);
static Template(xyrange_fn);
static Template(sflag_fn);
static Template(selfant_fn);
static Template(hist_fn);
static Template(uvsel_fn);
static Template(wtscal_fn);
static Template(peak_fn);
static Template(pwin_fn);
static Template(modfit_fn);
static Template(edmod_fn);
static Template(cpplt_fn);
static Template(addhis_fn);
static Template(clrhis_fn);
static Template(scangap_fn);
static Template(munit_fn);
static Template(addmc_fn);
static Template(uvstat_fn);
static Template(imstat_fn);
static Template(setcont_fn);
static Template(mapcol_fn);
static Template(mapfun_fn);
static Template(showpar_fn);
static Template(beamloc_fn);
static Template(polvec_fn);
static Template(specpl_fn);
static Template(specb_fn);
static Template(specp_fn);
static Template(spect_fn);
static Template(specuv_fn);
static Template(specop_fn);
static Template(specsm_fn);
static Template(specso_fn);
static Template(mapval_fn);
static Template(shiftto_fn);
static Template(read_models_fn);
static Template(write_models_fn);
static Template(clear_models_fn);
static Template(multi_model_fn);
static Template(mark_radec_fn);
static Template(mark_xy_fn);
static Template(clear_markers_fn);
static Template(wmarkers_fn);
static Template(delmarker_fn);
static Template(vis_stats_fn);
static Template(planet_temp_fn);
static Template(planet_geom_fn);
static Template(mjd_fn);
static Template(antenna_beam_fn);
static Template(pointing_center_fn);
static Template(primary_beam_fn);
static Template(flag_fn);
static Template(unflag_fn);
static Template(map_to_rad_fn);
static Template(rad_to_map_fn);
static Template(uv_to_wav_fn);
static Template(wav_to_uv_fn);

/*
 * Declare the function types below.
 */

static Functype dmapf_type[] = {
   {newob_fn,   NORM, 1,3,   " Cfl",    " 000",    " vvv",  1 },
   {mapsize_fn, NORM, 0,4,   " ifif",   " 0000",   " vvvv", 1 },
   {invert_fn,  NORM, 0,0,   "  ",      "  ",      "  ",    1 },
   {uvtaper_fn, NORM, 0,2,   " ff",     " 00",     " vv",   1 },
   {uvrange_fn, NORM, 0,2,   " ff",     " 00",     " vv",   1 },
   {uvwgt_fn,   NORM, 0,3,   " ffl",    " 000",    " vvv",  1 },
   {uvzero_fn,  NORM, 0,2,   " ff",     " 00",     " vv",   1 },
   {maplot_fn,  NORM, 0,2,   " Cl",     " 00",     " vv",   1 },
   {maplot_fn,  NORM, 0,2,   " Cl",     " 00",     " vv",   1 },
   {clean_fn,   NORM, 0,3,   " iff",    " 000",    " vvv",  1 },
   {restore_fn, NORM, 0,5,   " fffll",  " 00000",  " vvvvv",1 },
   {wbeam_fn,   NORM, 1,1,   " C",      " 0",      " v",    1 },
   {wmap_fn,    NORM, 1,1,   " C",      " 0",      " v",    1 },
   {wdmap_fn,   NORM, 1,1,   " C",      " 0",      " v",    1 },
   {wobs_fn,    NORM, 1,2,   " Cl",     " 00",     " vv",   1 },
   {wmodel_fn,  NORM, 0,2,   " Cl",     " 00",     " vv",   1 },
   {rmodel_fn,  NORM, 1,2,   " Cl",     " 00",     " vv",   1 },
   {gscal_fn,   NORM, 0,1,   " l",      " 0",      " v",    1 },
   {keep_fn,    NORM, 0,0,   "  ",      "  ",      "  ",    1 },
   {clrmod_fn,  NORM, 0,3,   " lll",    " 000",    " vvv",  1 },
   {shift_fn,   NORM, 0,2,   " ff",     " 00",     " vv",   1 },
   {unshift_fn, NORM, 0,0,   " ",       " ",       " ",     1 },
   {uvav_fn,    NORM, 1,2,   " fl",     " 00",     " vv",   1 },
   {head_fn,    NORM, 0,0,   "  ",      "  ",      "  ",    1 },
   {uncal_fn,   NORM, 0,3,   " lll",    " 000",    " vvv",  1 },
   {corpl_fn,   NORM, 0,2,   " Ci",     " 00",     " vv",   1 },
   {tname_fn,   NORM, 1,2,   "cii",     "000",     "vvv",   1 },
   {ntel_fn,    NORM, 0,1,   "ii",      "00",      "vv",    1 },
   {bname_fn,   NORM, 1,2,   "cii",     "000",     "vvv",   1 },
   {nbase_fn,   NORM, 0,1,   "ii",      "00",      "vv",    1 },
   {nsub_fn,    NORM, 0,0,   "i",       "0",       "v",     1 },
   {nif_fn,     NORM, 0,0,   "i",       "0",       "v",     1 },
   {nchan_fn,   NORM, 0,0,   "i",       "0",       "v",     1 },
   {addwin_fn,  NORM, 4,4,   " ffff",   " 0000",   " vvvv", 1 },
   {delwin_fn,  NORM, 0,0,   " ",       " ",       " ",     1 },
   {winmod_fn,  NORM, 0,1,   " l",      " 0",      " v",    1 },
   {startmod_fn,NORM, 0,2,   " Cf",     " 00",     " vv",   1 },
   {uvrad_fn,   NORM, 0,8,   " Cffffffl",  " 00000000",  " vvvvvvvv",  1 },
   {uvprj_fn,   NORM, 0,9,   " fCffffffl", " 000000000", " vvvvvvvvv", 1 },
   {self_fn,    NORM, 0,3,   " llf",    " 000",    " vvv",  1 },
   {staper_fn,  NORM, 0,2,   " ff",     " 00",     " vv",   1 },
   {slims_fn,   NORM, 0,2,   " ff",     " 00",     " vv",   1 },
   {rwins_fn,   NORM, 1,1,   " C",      " 0",      " v",    1 },
   {wwins_fn,   NORM, 0,2,   " Cl",     " 00",     " vv",   1 },
   {vplot_fn,   NORM, 0,4,   " iCii",   " 0000",   " vvvv", 1 },
   {uvplt_fn,   NORM, 0,4,   " Cffl",   " 0000",   " vvvv", 1 },
   {timpl_fn,   NORM, 0,2,   " Ci",     " 00",     " vv",   1 },
   {resof_fn,   NORM, 0,1,   " C",      " 0",      " v",    1 },
   {unoff_fn,   NORM, 0,2,   " ll",     " 00",     " vv",   1 },
   {save_fn,    NORM, 1,1,   " C",      " 0",      " v",    1 },
   {get_fn,     NORM, 1,1,   " C",      " 0",      " v",    1 },
   {loglev_fn,  NORM, 1,3,   " fff",    " 000",    " vvv",  1 },
   {xyrange_fn, NORM, 0,4,   " ffff",   " 0000",   " vvvv", 1 },
   {sflag_fn,   NORM, 0,3,   " lii",    " 000",    " vvv",  1 },
   {selfant_fn, NORM, 0,3,   " Clf",    " 000",    " vvv",  1 },
   {hist_fn,    NORM, 0,0,   " ",       " ",       " ",     1 },
   {uvsel_fn,   NORM, 0,MAXARG,        " Ci",     " 00",   " vv", 1 },
   {wtscal_fn,  NORM, 0,1,   "ff",      "00",      "?v",    1 },
   {peak_fn,    NORM, 1,2,   "fCC",     "000",     "vvv",   1 }, 
   {pwin_fn,    NORM, 0,2,   " fl",     " 00",     " vv",   1 },
   {modfit_fn,  NORM, 1,1,   " i",      " 0",      " v",    1 },
   {edmod_fn,   NORM, 0,1,   " l",      " 0",      " v",    1 },
   {cpplt_fn,   NORM, 0,4,   " iCii",   " 0000",   " vvvv", 1 },
   {addhis_fn,  NORM, 1,1,   " c",      " 0",      " v",    1 },
   {clrhis_fn,  NORM, 0,0,   " ",       " ",       " ",     1 },
   {scangap_fn, NORM, 0,2,   " fi",     " 00",     " vv",   1 },
   {munit_fn,   NORM, 0,1,   " C",      " 0",      " v",    1 },
   {addmc_fn,   NORM,4,15,   " flfflflflfliffl"," 000000000000000"," vvvvvvvvvvvvvvv",1},
   {uvstat_fn,  NORM, 1,1,   "fC",      "00",      "vv",    1 },
   {imstat_fn,  NORM, 1,1,   "fC",      "00",      "vv",    1 },
   {setcont_fn, NORM, 0,0,   " ",       " ",       " ",     1 },
   {mapcol_fn,  NORM, 0,3,   " Cff",    " 000",    " vvv",  1 }, 
   {mapfun_fn,  NORM, 0,3,   " Cff",    " 000",    " vvv",  1 },
   {showpar_fn, NORM, 0,0,   " ",       " ",       " ",     1 },
   {beamloc_fn, NORM, 0,4,   " ffff",   " 0000",   " vvvv", 1 },
   {polvec_fn,  NORM, 0,5,   " fffii",  " 00000",  " vvvvv",1 },
   {specpl_fn,  NORM, 0,9,   " iCiiffffi"," 000000000"," vvvvvvvvv", 1 },
   {specb_fn,   NORM, 0,MAXARG,  " C",  " 0",      " v",    1 },
   {specp_fn,   NORM, 0,MAXARG,  " C",  " 0",      " v",    1 },
   {spect_fn,   NORM, 0,3,   " CCf",    " 000",    " vvv",  1 },
   {specuv_fn,  NORM, 0,3,   " fff",    " 000",    " vvv",  1 },
   {specop_fn,  NORM, 0,2,   " CC",     " 00",     " vv",   1 },
   {specsm_fn,  NORM, 0,3,   " CCf",    " 000",    " vvv",  1 },
   {specso_fn,  NORM, 0,3,   " CCC",    " 000",    " vvv",  1 },
   {mapval_fn,  NORM, 2,2,   "fff",     "000",     "vvv",   1 },
   {shiftto_fn, NORM, 2,2,   " CC",     " 00",     " vv",   1 },
   {read_models_fn,  NORM, 1,1,  " C",  " 0",      " v",    1 },
   {write_models_fn, NORM, 1,1,  " C",  " 0",      " v",    1 },
   {clear_models_fn, NORM, 0,0,  " ",   " ",       " ",     1 },
   {multi_model_fn,  NORM, 0,1,  " l",  " 0",      " v",    1 },
   {mark_radec_fn,   NORM, 3,9,  " CCCifcfff", " 000000000", " vvvvvvvvv", 1 },
   {mark_xy_fn,      NORM, 3,9,  " ffCifcfff", " 000000000", " vvvvvvvvv", 1 },
   {clear_markers_fn,NORM, 0,0,  " ",   " ",       " ",     1 },
   {wmarkers_fn,     NORM, 0,1,  " C",  " 0",      " v",    1 },
   {delmarker_fn,    NORM, 2,2,  " CC", " 00",     " vv",   1 },
   {vis_stats_fn,    NORM, 1,1,  "fC",  "10",      "?v",    1 },
   {planet_temp_fn,  NORM, 1,4,  "ffCff", "00000", "?vvvv", 1 },
   {planet_geom_fn,  NORM, 1,2,  "fCf",  "100",    "?vv",   1 },
   {mjd_fn,          NORM, 1,1,  "fC",   "00",     "vv",    0 },
   {antenna_beam_fn, NORM, 1,4,  " Cfff"," 0100",  " vvvv", 1 },
   {pointing_center_fn, NORM, 0,2, " CC"," 00",    " vv",   1 },
   {primary_beam_fn, NORM, 0,3,  " fff", " 100",   " vvv",  1 },
   {flag_fn,         NORM, 1,4,  " ClCC"," 0000",  " vvvv", 1 },
   {unflag_fn,       NORM, 1,4,  " ClCC"," 0000",  " vvvv", 1 },
   {map_to_rad_fn,   NORM, 1,1,  "ff",   "00",     "vv",    0 },
   {rad_to_map_fn,   NORM, 1,1,  "ff",   "00",     "vv",    0 },
   {uv_to_wav_fn,    NORM, 1,1,  "ff",   "00",     "vv",    0 },
   {wav_to_uv_fn,    NORM, 1,1,  "ff",   "00",     "vv",    0 },
};

/*
 * In the same order as the above array of types, define the array
 * of user names for the functions.
 */

static char *dmapf_name[] = {
   "observe",
   "mapsize",
   "invert",
   "uvtaper",
   "uvrange",
   "uvweight",
   "uvzero",
   "mapplot",
   "maplot",
   "clean",
   "restore",
   "wbeam",
   "wmap",
   "wdmap",
   "wobs",
   "wmodel",
   "rmodel",
   "gscale",
   "keep",
   "clrmod",
   "shift",
   "unshift",
   "uvaver",
   "header",
   "uncalib",
   "corplot",
   "telname",
   "ntel",
   "basename",
   "nbase",
   "nsub",
   "nif",
   "nchan",
   "addwin",
   "delwin",
   "winmod",
   "startmod",
   "radplot",
   "projplot",
   "selfcal",
   "selftaper",
   "selflims",
   "rwins",
   "wwins",
   "vplot",
   "uvplot",
   "tplot",
   "resoff",
   "clroff",
   "save",
   "get",
   "loglevs",
   "xyrange",
   "selfflag",
   "selfant",
   "showhist",
   "select",
   "wtscale",
   "peak",
   "peakwin",
   "modelfit",
   "edmodel",
   "cpplot",
   "addhist",
   "clrhist",
   "scangap",
   "mapunits",
   "addcmp",
   "uvstat",
   "imstat",
   "setcont",
   "mapcolor",
   "mapfunc",
   "showpar",
   "beamloc",
   "polvec",
   "specplot",
   "specbase",
   "specpol",
   "spectime",
   "specuvr",
   "specopt",
   "specsmooth",
   "specorder",
   "mapvalue",
   "shiftto",
   "read_models",
   "write_models",
   "clear_models",
   "multi_model",
   "mark_radec",
   "mark_xy",
   "clear_markers",
   "wmarkers",
   "delmarker",
   "vis_stats",
   "planet_temp",
   "planet_geometry",
   "mjd",
   "antenna_beam",
   "pointing_center",
   "primary_beam",
   "flag",
   "unflag",
   "map_to_rad",
   "rad_to_map",
   "uv_to_wav",
   "wav_to_uv",
};

/*
 * List special help topics.
 */
static char *dmap_help[] = {
  "whatsnew",
  "spectral_line",
  "models",
  "editing",
  "subarrays",
  "multi_if",
  "polarization",
  "antenna_names",
};

static int dmap_begin(void);  /* Module initialization function */
static EXITFN(dmap_end);      /* Module closedown function */

/*
 * Record the above declarations etc for this module in a global
 * structure for use when building the main symbol table.
 */

Module m_difmap = {
  "difmap",
  HELP_DIR,
  dmap_help, COUNT(dmap_help),
  dmapv_type, dmapv_name, COUNT(dmapv_name),
  dmapf_type, dmapf_name, COUNT(dmapf_name),
  dmap_begin,
  dmap_end
};

/*
 * External signal detection variable.
 */
extern int no_error;

/*
 * Local variables and functions.
 */

static Observation *vlbob = NULL;
static MapBeam *vlbmap    = NULL;
static Mapwin *vlbwins    = NULL;
static Specattr *vlbspec  = NULL;
static MarkerList *mapmarkers = NULL;

/*
 * Assign suffixes to file types.
 */
enum{MAXSUF=6};    /* The max suffix length including \0 */
static const char *uvf_nam=".uvf";   /* Merge file suffix */
static const char *mod_nam=".mod";   /* Model file suffix */
static const char *cmod_nam=".cmod"; /* Continuum-model file suffix */
static const char *win_nam=".win";   /* Window file suffix */
static const char *fits_nam=".fits"; /* Window file suffix */
static const char *par_nam=".par";   /* Parameter command file suffix */
static const char *mtab_nam=".mtab"; /* A model-table file */

static void obs_end();
static int wrtpars(char *parname, char *basename);
static int w_flt_array(FILE *fp, char *name, Descriptor *dsc);
static int write_marker_commands(FILE *fp);

static int nodata(const char *fname, Obstate state);
static int nomap(const char *fname);

/*
 * Enumerate the values of vlbmap->domap.
 */
enum {
  MAP_IS_MAP=0, /* The map array contains the latest residual map */
  MAP_IS_STALE, /* The map array doesn't contain anything useful */
  MAP_IS_CLEAN, /* The map array contains the latest clean map */
  MAP_IS_PMAP,  /* The map array contains the latest residual map, along the */
                /*  latest residual polarization intensity and angle maps */
                /*  overwriting the margins above and below the main map. */
  MAP_IS_PCLN   /* The map array contains the latest restored map plus the */
                /*  latest restored polarization intensity and angle maps */
                /*  overwriting the margins above and below the main map. */
} MapState;

/*.......................................................................
 * Difmap module specific startup code.
 *
 * Output:
 *  return    int   0 - OK.
 *                 -1 - Error.
 */
static int dmap_begin(void)
{
/*
 * Start a log file.
 */
  logfile("difmap.log");
/*
 * Set up the mapplot color table.
 */
  mappar.ctab = new_Ctable();
  if(mappar.ctab == NULL)
    return -1;
/*
 * Create an empty table of map markers.
 */
  mapmarkers = new_MarkerList();
  if(!mapmarkers)
    return -1;
  return 0;
}

/*.......................................................................
 * Difmap module specific closedown code.
 *
 * Input:
 *  code    Exitcode   DO_QUIT  -  Minimal fast closedown.
 *                     DO_EXIT  -  Full closedown.
 */
static EXITFN(dmap_end)
{
  int waserr=0;       /* Error status while saving files */
  enum {MAX_TRY=5};   /* Max number of attempts to get file name prefix */
  int try_count=0;    /* Number of attempts to get file name prefix */
  enum {MAX_PRE=80};  /* Max length of filename prefix string */
  char *cptr;         /* Pointer into reply[] */
  static char reply[MAX_PRE-1]; /* User's reply */
  static Descriptor filearg = {'c',0,NO_DEL,1,{1,1,1},NULL};
  static Descriptor *invals = &filearg;  /* Argument to send to save_fn() */
  char *prompt;       /* The prompt to use to query the user */
/*
 * Close PGPLOT.
 */
  cpgend();
/*
 * If full closedown has been requested give the user the option
 * of having their files saved.
 */
  if(code==DO_EXIT && vlbob!=NULL) {
/*
 * Prompt the user for a file name prefix to save the files to
 * until they enter a valid prefix, until they press return without
 * enterring a prefix, or the trial count expires.
 */
    do {
/*
 * Prompt the user for a file name prefix.
 */
      switch(try_count) {
      case 0:
	prompt =
	 "Enter a file name prefix, or press return to quit without saving: ";
	break;
      case MAX_TRY-1:
	prompt = "This is your last chance to enter a prefix: ";
	break;
      default:
	prompt = "Try a different prefix: ";
	break;
      };
/*
 * Read the prefix.
 */
      waserr = lexgets(reply, MAX_PRE, stdin, prompt) != 0;
/*
 * Skip white-space to see if a prefix was provided.
 */
      if(!waserr) {
	for(cptr=reply; isspace((int)*cptr); cptr++);
	if(*cptr != '\0') {
	  VOIDPTR(&filearg) = &cptr;
	  waserr = save_fn(&invals, 1, NULL);
	};
      };
    } while(waserr && ++try_count<MAX_TRY);
  };
/*
 * Clean up memory and scratch files.
 */
  obs_end();
/*
 * Delete colormap symbol table.
 */
  mappar.ctab = del_Ctable(mappar.ctab);
/*
 * Delete the map and beam.
 */
  vlbmap = del_MapBeam(vlbmap);
/*
 * Delete the list of map markers.
 */
  mapmarkers = del_MarkerList(mapmarkers);
/*
 * Close the log file.
 */
  logfile(NULL);
  return;
}

/*.......................................................................
 * Read in a new vlbi observation.
 *
 * Input:
 *  name      char *  The name of the UV FITS file.
 *  binwid   float    The integration time to the bin the visibilities
 *                    to (seconds). No binning will be performed if
 *                    binwid<1.0 second.
 *  scatter   char    If TRUE replace the data weights by weights
 *                    derived from the data scatter.
 */
static Template(newob_fn)
{
  char hisline[81];  /* String to compose history info in */
  char *cptr;        /* Pointer to return value of date_str() */
  char *name="";     /* The file name of the requested merge file */
  double binwid=0.0; /* The integration bin width (seconds) */
  int scatter=0;     /* True if weights are to be found from the data scatter */
/*
 * Get the arguments.
 */
  switch(npar) {
  case 3:
    scatter = *LOGPTR(invals[2]);
  case 2:
    binwid = *FLTPTR(invals[1]);
  case 1:
    name = *STRPTR(invals[0]);
  };
/*
 * Check that the requested merge file exists before initializing.
 */
  if(!file_exists(name)) {
    lprintf(stderr, "observe: File \"%s\" does not exist\n", name);
    return -1;
  };
/*
 * If we already have an observation loaded, delete it and its scratch
 * files, and any associated data structures including models.
 */
  obs_end();
/*
 * Clear the map and beam.
 */
  if(vlbmap)
    vlbmap = new_MapBeam(vlbmap, vlbmap->nx, vlbmap->xinc,
			 vlbmap->ny, vlbmap->yinc);
/*
 * Discard the current list of markers.
 */
  clr_MarkerList(mapmarkers);
/*
 * Read the new merge file into an Observation struct.
 */
  vlbob = new_Observation(name, binwid, scatter, 1, NULL, NO_POL);
  if(vlbob == NULL)
    return -1;
/*
 * Allocate a new specplot attributes descriptor.
 */
  vlbspec = new_Specattr(vlbob);
  if(!vlbspec) {
    obs_end();
    return -1;
  };
/*
 * Reset selected parameters.
 */
  invpar = invdef;
  respar = resdef;
  slfpar = slfdef;
/*
 * Append a history line so that if this structure is written to a
 * new FITS file there will be some evidence that difmap got its
 * grubby hands on it.
 */
  cptr = date_str();
  if(cptr != NULL) {
    sprintf(hisline, "DIFMAP  Read into difmap on %.*s", 48, cptr);
    add_hist(vlbob, hisline);
  };
/*
 * File read OK.
 */
  return no_error;
}

/*.......................................................................
 * Free up the dynamically allocated memory associated with the current
 * observation.
 */
static void obs_end()
{
/*
 * Delete the observation descriptor.
 */
  if(vlbob)
    vlbob = del_Observation(vlbob);
/*
 * Delete the list of CLEAN windows.
 */
  vlbwins = del_Mapwin(vlbwins);
/*
 * Delete the specplot attributes descriptor of this observation.
 */
  vlbspec = del_Specattr(vlbspec);
  return;
}

/*.......................................................................
 * Create a new MapBeam instance for a given mapsize and cellsize.
 *
 * Input:
 *  nx     int  Number of pixels along X-axis of map/beam grid (power of 2).
 *  xinc float  The X-axis cellsize of the grid in map xy units.
 *  ny     int  Number of pixels along X-axis of map/beam grid (power of 2).
 *  yinc float  The Y-axis cellsize of the grid in map xy units.
 */
static Template(mapsize_fn)
{
/*
 * Are we changing the map size?
 */
  if(npar > 0) {
    float xinc;    /* Copy of passed X-axis cellsize converted to radians */
    float yinc;    /* Copy of passed Y-axis cellsize converted to radians */
    int nx;        /* Copy of passed X-axis dimension */
    int ny;        /* Copy of passed Y-axis dimension */
/*
 * Two few parameters?
 */
    if(npar<1) {
      lprintf(stderr, "mapsize: X-axis number of pixels required.\n");
      return -1;
    };
/*
 * Get the user arguments. The X-axis grid size and cell-size are mandatory.
 */
    nx = *INTPTR(invals[0]);
    xinc = npar>1 ? xytorad(fabs(*FLTPTR(invals[1]))) : 0.0;
/*
 * The second two are optional - substitute the X-axis values if none
 * are given.
 */
    ny =   npar>2 ? *INTPTR(invals[2]) : nx;
    yinc = npar>3 ? xytorad(fabs(*FLTPTR(invals[3]))) : xinc;
/*
 * Fill in the optimal pixel sizes if no pixel size was given.
 */
    if(xinc <= 0.0 || yinc <= 0.0) {
      float xmax, ymax;
      if(nodata("mapsize", OB_SELECT) ||
	 optimal_pixel_size(vlbob, invpar.uvmin, invpar.uvmax, nx, ny,
			    &xmax, &ymax))
	return -1;
/*
 * If a square map is needed, substitute the minimum of the maximum
 * x and y pixel sizes.
 */
      if(npar <= 3) {
	xinc = yinc = xmax < ymax ? xmax : ymax;
      } else {
	if(xinc <= 0.0)
	  xinc = xmax;
	if(yinc <= 0.0)
	  yinc = ymax;
      };
    };
/*
 * Get a new map and beam (reusing the old one if possible).
 */
    vlbmap = new_MapBeam(vlbmap, nx, xinc, ny, yinc);
    if(vlbmap == NULL)
      return -1;
/*
 * Make the user 'map' and 'beam' variables point at the new map and
 * beam.
 */
    VOIDPTR(&mb_beam) = vlbmap->beam;
    VOIDPTR(&mb_map) = vlbmap->map;
    mb_beam.adim[0] = nx;
    mb_beam.adim[1] = ny;
    mb_map.adim[0] = nx;
    mb_map.adim[1] = ny;
    mb_beam.num_el = mb_map.num_el = (size_t) nx*ny;
  };
/*
 * Report the curren map dimensions.
 */
  if(vlbmap==NULL) {
    lprintf(stdout, "No map has been allocated yet.\n");
  } else {
    lprintf(stdout,
	    "Map grid = %dx%d pixels with %#.3gx%#.3g %s cellsize.\n",
	    vlbmap->nx, vlbmap->ny, radtoxy(vlbmap->xinc),
	    radtoxy(vlbmap->yinc), mapunits(U_TLAB));
  };
  return no_error;
}

/*.......................................................................
 * Specify a UV gaussian taper for weighting down long baselines.
 * If the specified gaussian value is not between 0 and 1 or the
 * specified radius is <= 0 then cancel any previous taper.
 *
 * Input:
 *  gauval    float  The gaussian taper value at 'gaurad'.
 *  gaurad    float  The radius (user UV distance units) that the taper
 *                   has value 'gauval'.
 */
static Template(uvtaper_fn)
{
/*
 * The UV taper is undefined until an observation has been read.
 */
  if(nodata("uvtaper", OB_INDEX))
    return -1;
/*
 * Is anything being changed?
 */
  if(npar > 0) {
/*
 * Override current taper values.
 */
    switch(npar) {
    case 2:
      invpar.gaurad = uvtowav(*FLTPTR(invals[1]));  /* Convert to wavelengths */
    case 1:
      invpar.gauval = *FLTPTR(invals[0]);
    };
/*
 * The change in weighting will require the map and beam to be recomputed.
 */
    if(vlbmap)
      vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
/*
 * Ensure that the new estimated clean-beam size is substituted on the
 * next restore.
 */
    respar.doauto = 1;
/*
 * Cancel the taper if the new values are out of bounds.
 */
    if(invpar.gauval<=0.0f || invpar.gauval>=0.99f || invpar.gaurad<=0.0f) {
      invpar.gauval = 0.0f;
      invpar.gaurad = 0.0f;
    };
  };
/*
 * Report the current settings.
 */
  if(invpar.gauval <= 0.0f || invpar.gaurad <= 0.0f) {
    lprintf(stdout, "No UV-taper is currently set.\n");
  } else {
    lprintf(stdout,
	   "Gaussian taper: value %g at UV radius = %g %s.\n",
	   invpar.gauval, wavtouv(invpar.gaurad), uvwunits(U_TLAB));
  };
  return no_error;
}

/*.......................................................................
 * Set or unset a UV cutoff.
 *
 * Input:
 *   uvmin  float   The min UV radius cutoff in user UV distance units.
 *   uvmax  float   The max UV radius cutoff in user UV distance units.
 *                  If the largest of uvmin and uvmax <= 0.0, no cuttoff
 *                  will be selected.
 */
static Template(uvrange_fn)
{
/*
 * The UV range is undefined until an observation has been read.
 */
  if(nodata("uvrange", OB_INDEX))
    return -1;
/*
 * Is anything being changed?
 */
  if(npar > 0) {
/*
 * Override defaults.
 */
    switch(npar) {
    case 2:
      invpar.uvmax = uvtowav(*FLTPTR(invals[1]));
    case 1:
      invpar.uvmin = uvtowav(*FLTPTR(invals[0]));
    };
/*
 * The change in uvrange will require the map and beam to be recomputed.
 */
    if(vlbmap)
      vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
/*
 * Ensure that the new estimated clean-beam size is substituted on the
 * next restore.
 */
    respar.doauto = 1;
/*
 * Cancel the UV range if any values are out of bounds.
 */
    if(invpar.uvmin >= invpar.uvmax || invpar.uvmax <= 0.0f)
      invpar.uvmin = invpar.uvmax = 0.0f;
/*
 * Enforce positivity on uvmin.
 */
    if(invpar.uvmin < 0.0f)
      invpar.uvmin = 0.0f;
  };
/*
 * Inform user of the interpretation of their selection.
 */
  if(invpar.uvmax > 0.0f)
    lprintf(stdout,
	    "Only data in the UV range: %g -> %g (%s) will be gridded.\n",
	    wavtouv(invpar.uvmin), wavtouv(invpar.uvmax), uvwunits(U_TLAB));
  else
    lprintf(stdout, "The full UV range of the data is currently selected for gridding.\n");
  return no_error;
}

/*.......................................................................
 * Set (or unset) a zero baseline flux.
 *
 * Input:
 *  uvzero  float   The zero-baseline flux (Jy).
 *  weight  float   The visibility weight to assign the zero-baseline
 *                  flux.
 */
static Template(uvzero_fn)
{
  float weight;  /* Zero-baseline flux visibility weight */
  float flux;    /* Zero-baseline flux */
/*
 * The zero-spacing flux is undefined until an observation has been read.
 */
  if(nodata("uvzero", OB_INDEX))
    return -1;
/*
 * Is the flux to be changed?
 */
  if(npar > 0) {
/*
 * Get the visibility weight to assign to the zero-baseline flux.
 */
    if(npar >= 2) {
      weight = *FLTPTR(invals[1]);
    } else if(vlbob->uvzero.wt > 0.0f) {
      weight = vlbob->uvzero.wt;
    } else {
      lprintf(stderr,
	 "uvzero: Warning - substituting 1.0 for missing visibility weight.\n");
      weight = 1.0f;
    };
/*
 * Get the zero-baseline flux.
 */
    flux = *FLTPTR(invals[0]);
/*
 * Record the new values.
 */
    if(weight > 0.0f) {
      vlbob->uvzero.amp = flux;
      vlbob->uvzero.wt = weight;
    } else {
      vlbob->uvzero.amp = 0.0f;
      vlbob->uvzero.wt = 0.0f;
    };
/*
 * The current map and beam are now invalid.
 */
    if(vlbmap)
      vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
/*
 * Ensure that the new estimated clean-beam size is substituted on the
 * next restore.
 */
    respar.doauto = 1;
  };
/*
 * Inform user of the interpretation of their selection.
 */
  if(vlbob->uvzero.wt > 0.0f)
    lprintf(stdout, "Zero-baseline flux set to %g Jy. Weight=%g\n",
	    vlbob->uvzero.amp, vlbob->uvzero.wt);
  else
    lprintf(stdout, "Zero-baseline flux not set.\n");
  return no_error;
}

/*.......................................................................
 * Set up the UV-grid weighting options (apart from tapering).
 *
 * Input:
 *   uvbin  float   The bin size for uniform weighting in pixels. If <=0
 *                  uniform weighting is cancelled.
 *   errpow float   If < 0.0 then the amplitude errors, raised to the
 *                  power 'errpow', will be used to scale the weights.
 *   dorad  char    If TRUE then radial weighting is to be applied.
 */
static Template(uvwgt_fn)
{
/*
 * The weight is undefined until an observation has been read.
 */
  if(nodata("uvweight", OB_INDEX))
    return -1;
/*
 * Is anything to be changed?
 */
  if(npar > 0) {
/*
 * Set the weighting according to the number of arguments sent.
 */
    switch (npar) {
    case 3:
      invpar.dorad = *LOGPTR(invals[2]);
    case 2:                               /* Note fallthrough from case 3: */
      invpar.errpow = *FLTPTR(invals[1]);
      if(invpar.errpow > 0.0f)
	invpar.errpow = 0.0f;
    case 1:                               /* Note fallthrough from case 2: */
      invpar.uvbin  = *FLTPTR(invals[0]);
      if(invpar.uvbin < 0.0f)
	invpar.uvbin = 0.0f;
      else if(invpar.uvbin > 0.0f && invpar.uvbin < 1.0)
	invpar.uvbin = 1.0f;
    };
/*
 * The change in weighting invalidates the current map and beam.
 */
    if(vlbmap)
      vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
/*
 * Ensure that the new estimated clean-beam size is substituted on the
 * next restore.
 */
    respar.doauto = 1;
  };
/*
 * Report the weighting selection.
 */
  if(invpar.uvbin > 0.0f) {
    lprintf(stdout, "Uniform weighting binwidth: %g (pixels).\n",invpar.uvbin);
  } else {
    lprintf(stdout, "Uniform weighting is not currently selected.\n");
  };
  if(invpar.errpow < 0.0f) {
    lprintf(stdout,
        "Gridding weights will be scaled by errors raised to the power %g.\n",
	invpar.errpow);
  } else {
    lprintf(stdout, "Amplitude error weighting is not currently selected.\n");
  };
  lprintf(stdout, "Radial weighting is %scurrently selected.\n",
	  invpar.dorad ? "" : "not ");
  return no_error;
}

/*.......................................................................
 * Grid the UV data currently in the 'vlbob' observation structure
 * and fourier invert it to generate a dirty map and dirty beam.
 */
static Template(invert_fn)
{
/*
 * Check that a map and beam have been allocated and that an observation
 * has been read.
 */
  if(nomap("invert") || nodata("invert", OB_SELECT))
    return -1;
/*
 * Invert the data.
 */
  if(uvinvert(vlbob, vlbmap, invpar.uvmin, invpar.uvmax, invpar.gauval,
	      invpar.gaurad, invpar.dorad, invpar.errpow,
	      invpar.uvbin))
    return -1;
/*
 * Record the latest estimate of the equivalent restoring-beam size, for
 * use by the next restore.
 */
  respar.e_bmin = vlbmap->e_bmin;
  respar.e_bmaj = vlbmap->e_bmaj;
  respar.e_bpa  = vlbmap->e_bpa  * rtod;
  return no_error;
}

/*.......................................................................
 * Clean a dirty map.
 *
 * Input:
 *   niter    int    The max number of iterations.
 *   gain   float    The CLEAN loop gain < 1.0 .
 *   cutoff float    The residual flux to stop cleaning at. (Jy/Beam).
 */
static Template(clean_fn)
{
  Model *clnmod;   /* The clean model */
/*
 * Check that a map and beam have been allocated and that we have
 * something to map.
 */
  if(nodata("clean", OB_SELECT) || nomap("clean"))
    return -1;
/*
 * Re-invert the data if necessary.
 */
  if((vlbmap->domap || vlbmap->dobeam) && invert_fn(invals,0,outvals))
    return -1;
/*
 * Over-ride default clean parameters where arguments are given.
 * Note fallthrough between cases (no 'break' statements).
 */
  switch(npar) {
  case 3:
    clnpar.cutoff = *FLTPTR(invals[2]);
  case 2:
    clnpar.gain = *FLTPTR(invals[1]);
  case 1:
    clnpar.niter = *INTPTR(invals[0]);
  };
/*
 * Report current operating paremeter settings.
 */
  lprintf(stdout, "clean: niter=%d  gain=%g  cutoff=%g\n",
	  clnpar.niter, clnpar.gain, clnpar.cutoff);
/*
 * CLEAN the map.
 */
  clnmod = mapclean(vlbob, vlbmap, vlbwins, clnpar.niter, clnpar.cutoff,
		    clnpar.gain, 1);
  if(!clnmod)
    return -1;
/*
 * Correct for the effects of the combined primary beams?
 */
  if(count_antenna_beams(vlbob->ab) > 0) {
    Modcmp *cmp;
    for(cmp=clnmod->head; cmp; cmp=cmp->next)
      pb_correct_delta_cmp(vlbob, cmp);
  };
/*
 * Merge the new clean model with the existing model.
 */
  add_mod(vlbob->newmod, clnmod, 1, 1);
/*
 * Report the combined flux in the two models.
 */
  lprintf(stdout, "Combined flux in latest and established models = %g Jy\n",
	 vlbob->newmod->flux + vlbob->model->flux);
  return no_error;
}

/*.......................................................................
 * Restore a residual map using either the past model, or the 
 * sub-model generated in the last invokation of clean.
 *
 * Input:
 *  bmin  float  The length of the restoring beam's minor axis.
 *  bmaj  float  The length of the restoring beam's major axis.
 *  bpa   float  The restoring beam position angle of the major axis
 *               (North to East).
 *  noresid int  If true, clear the residual map before restoring the
 *               model.
 *  dosm    int  If true (default) smooth residuals before restoring.
 */
static Template(restore_fn)
{
  int noresid=0;  /* Default to restoring on top of residuals */
  int dosm=1;     /* Default to smoothing the residual map before restoring */
/*
 * Check that a map and beam have been allocated and that we have
 * something to map.
 */
  if(nodata("restore", OB_SELECT) || nomap("restore"))
    return -1;
/*
 * Re-invert the data if necessary, before restoring.
 */
  if(vlbmap->domap && invert_fn(invals,0,outvals))
    return -1;
/*
 * Override the default bmin,bmaj,bpa with any arguments that were sent.
 */
  switch(npar) {     /* NB. the case fallthroughs are intentional! */
  case 5:
    dosm = *LOGPTR(invals[4]);
  case 4:
    noresid = *LOGPTR(invals[3]);
  case 3:
    respar.bpa = *FLTPTR(invals[2]);
  case 2:
    respar.bmaj = xytorad(*FLTPTR(invals[1]));
  case 1:
    respar.bmin = xytorad(*FLTPTR(invals[0]));
  };
/*
 * Insert defaults for unspecified arguments.
 */
  switch(npar) { /* Case fallthrough deliberate */
  case 1:
    respar.bmaj = respar.bmin; /* Circular beam */
  case 2:
    respar.bpa = 0.0f;         /* Assume zero inclination angle */
  };
/*
 * If arguments were given, see if they should be used to over-ride the
 * automatic beam size determinations in this and future uses of restore.
 */
  respar.doauto = respar.bmin==0.0f || (respar.doauto && npar==0);
/*
 * Substitute the estimate from the last 'invert' if required.
 */
  if(respar.doauto) {
    lprintf(stdout,
      "restore: Substituting estimate of restoring beam from last 'invert'.\n");
    respar.bmin = respar.e_bmin;
    respar.bmaj = respar.e_bmaj;
    respar.bpa = respar.e_bpa;
  };
/*
 * Check the specified beam dimensions.
 */
  if(respar.bmin <= 0.0f) {
    lprintf(stderr,"restore: Illegal bmin=%g %s.\n",
	    radtoxy(respar.bmin), mapunits(U_TLAB));
    return -1;
  };
  if(respar.bmaj <= 0.0f) {
    lprintf(stderr,"restore: Illegal bmaj=%g %s.\n",
	    radtoxy(respar.bmaj), mapunits(U_TLAB));
    return -1;
  };
/*
 * Swap bmin and bmaj if bmin>bmaj.
 */
  if(respar.bmin>respar.bmaj) {
    float ftmp  = respar.bmin;
    respar.bmin = respar.bmaj;
    respar.bmaj = ftmp;
  };
/*
 * Report the beam size being used.
 */
  lprintf(stdout,
     "Restoring with beam: %.4g x %.4g at %.4g degrees (North through East)\n",
     radtoxy(respar.bmin), radtoxy(respar.bmaj), respar.bpa);
/*
 * Restore the map from the established model.
 */
  if(vlbob->model->ncmp + vlbob->newmod->ncmp < 1) {
    lprintf(stdout, "No model to restore with.\n");
    return -1;
  } else {
/*
 * Mark the map as invalid, in case of errors.
 */
    vlbmap->domap = MAP_IS_STALE;
#if 0
/*
 * First correct the residual map for the primary beam.
 */
    if(count_antenna_beams(vlbob->ab) > 0 && pb_cor_map(vlbob, vlbmap, 0.005))
      return 1;
#endif
/*
 * Restore established model.
 */
    if(vlbob->model->ncmp>0) {
      if(mapres(vlbob, vlbmap, vlbob->model, vlbmap->map, respar.bmaj,
		respar.bmin, respar.bpa*dtor, 0, noresid, dosm,
		getfreq(vlbob, -1)) == NULL)
	return -1;
/*
 * Prevent the following mapres() from smoothing the partially restored map.
 */
      dosm = 0;
    };
/*
 * Restore the tentative model.
 */
    if(vlbob->newmod->ncmp>0) {
      if(mapres(vlbob, vlbmap, vlbob->newmod, vlbmap->map, respar.bmaj,
		respar.bmin, respar.bpa*dtor, 0, noresid, dosm,
		getfreq(vlbob, -1)) == NULL)
	return -1;
    };
/*
 * Mark the residual map as invalid, while signalling that the restored
 * map is up to date.
 */
    vlbmap->domap = MAP_IS_CLEAN;
  };
  return no_error;
}

/*.......................................................................
 * Write a fits file of the restored map with a clean-component table
 * extension.
 *
 * Input:
 *   name       char *   The file name for the fits file.
 */
static Template(wmap_fn)
{
/*
 * Make sure that an observation has been read and that map and beam
 * have been allocated.
 */
  if(nodata("wmap", OB_SELECT) || nomap("wmap"))
    return -1;
/*
 * Only one model can be saved in the FITS file. Establish the
 * latest clean model.
 */
  if(keep_fn(invals,0,outvals) == -1)
    return -1;
/*
 * If the map is not already restored, restore it.
 */
  if(!vlbmap->ncmp || vlbmap->domap!=MAP_IS_CLEAN) {
    if(restore_fn(invals,0,outvals) == -1)
      return -1;
  };
/*
 * Write the file.
 */
  if(w_MapBeam(vlbob, vlbmap, 1, *STRPTR(invals[0])))
    return -1;
  return no_error;
}

/*.......................................................................
 * Write a fits file of the clean beam.
 *
 * Input:
 *   name       char *   The file name for the fits file.
 */
static Template(wbeam_fn)
{
/*
 * Make sure that an observation has been read and that map and beam
 * have been allocated.
 */
  if(nodata("wbeam", OB_SELECT) || nomap("wbeam"))
    return -1;
/*
 * If the beam is not up to date, re-invert.
 */
  if(vlbmap->dobeam && invert_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Write the file.
 */
  if(w_MapBeam(vlbob, vlbmap, 0, *STRPTR(invals[0])))
    return -1;
  return no_error;
}

/*.......................................................................
 * Write a fits file of the residual map.
 *
 * Input:
 *   name       char *   The file name for the fits file.
 */
static Template(wdmap_fn)
{
/*
 * Make sure that an observation has been read and that map and beam
 * have been allocated.
 */
  if(nodata("wdmap", OB_SELECT) || nomap("wdmap"))
    return -1;
/*
 * If the map is not up to date, re-invert.
 */
  if(vlbmap->domap && invert_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Write the file.
 */
  if(w_MapBeam(vlbob, vlbmap, 1, *STRPTR(invals[0])))
    return -1;
  return no_error;
}

/*.......................................................................
 * Write a UV FITS file of the current UV data.
 *
 * Input:
 *  name        char *  The name of the file to be written.
 *  doshift  logical    Shift the data in the output file according to
 *                      vlbob->geom.east and vlbob->geom.north, if doshift
 *                      is true.
 *                      Default=false.
 */
static Template(wobs_fn)
{
  char *filename = NULL;  /* The name to give the output file */
  int doshift = 0;        /* True to shift the pointing center in the */
                          /*  output file. */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("wobs", OB_INDEX))
    return -1;
/*
 * Get the arguments.
 */
  switch(npar) { /* Case fallthrough deliberate */
  case 2:
    doshift = *LOGPTR(invals[1]);
  case 1:
  default:
    filename = *STRPTR(invals[0]);
  };
/*
 * Write the file.
 */
  if(uvf_write(vlbob, filename, doshift))
    return -1;
  return no_error;
}

/*.......................................................................
 * Write the latest and established CLEAN models to a vlbi model file.
 *
 * Input:
 *  name   char *  The file name for the model file.
 *  docont bool    If true, write the continuum model only. If false or
 *                 omitted write the normal models only.
 */
static Template(wmodel_fn)
{
  char buf[80];   /* Buffer used while composing RA DEC header lines */
  int waserr=0;   /* Records whether wmodel() succeded */
  int docont=0;   /* True to write continuum */
  float east;     /* Any eastward offset to remove while writing */
  float north;    /* Any northward offset to remove while writing */
  Model *model;   /* Established model to be written */
  Model *newmod;  /* Tentative model to be written */
/*
 * Default to write to stdout.
 */
  FILE *modfd=stdout;     /* File descriptor for model file */
  char *modfil="(stdout)";/* Pointer to model-file name */
/*
 * Sanity check.
 */
  if(nodata("wmodel", OB_INDEX))
    return -1;
/*
 * Get arguments.
 */
  switch(npar) {
  case 2:
    docont = *LOGPTR(invals[1]);
  case 1:
/*
 * Override stdout default?
 */
    if(*STRPTR(invals[0])[0] != '\0') {
/*
 * Get the model-file name.
 */
      modfil = *STRPTR(invals[0]);
/*
 * Open the model file.
 */
      modfd = fopen(modfil, "w");
      if(modfd == NULL) {
	lprintf(stderr, "wmodel: Unable to open new model file: %s\n", modfil);
	return -1;
      };
    };
  };
/*
 * Get pointer to the models to be written.
 */
  model = docont ? vlbob->cmodel : vlbob->model;
  newmod = docont ? vlbob->cnewmod : vlbob->newmod;
/*
 * Report what is being done.
 */
  lprintf(stdout, "Writing %d %smodel components to file: %s\n",
	 (model->ncmp + newmod->ncmp), docont ? "continuum ":"", modfil);
/*
 * Determine the coordinate offsets to remove from the model components.
 */
  east  = vlbob->geom.east;
  north = vlbob->geom.north;
/*
 * Write the RA and DEC of the unshifted phase-center.
 */
  waserr = waserr ||
    lprintf(modfd,"! Center RA: %s,  ",
	    sradhms(vlbob->source.ra, 5, 0, buf))<0 ||
    lprintf(modfd, "Dec: %s (%.1f)\n", sraddms(vlbob->source.dec, 5, 0, buf),
	    vlbob->source.epoch) < 0;
/*
 * Write the established and tentative models.
 */
  if(model->ncmp > 0) {
    waserr = waserr || lprintf(modfd, "! Established model.\n") < 0;
    waserr = waserr || wmodel(model, east, north, 0, 0.0f, modfd);
  };
  if(newmod->ncmp > 0) {
    waserr = waserr || lprintf(modfd, "! Tentative model.\n") < 0;
    waserr = waserr || wmodel(newmod, east, north, 0, 0.0f, modfd);
  };
  if((modfd!=stdout && fclose(modfd)==EOF) || waserr) {
    lprintf(stderr, "wmodel: Error writing file: %s\n", modfil);
    return -1;
  };
  return no_error;
}

/*.......................................................................
 * Replace the latest clean model with the model in a specified VLBI
 * model file.
 *
 * Input:
 *  modfil  char *  The name of the model file.
 *  docont  bool    If true, read into the continuum model.
 */
static Template(rmodel_fn)
{
  char *modfil = "";  /* The name of the model file */
  float east;   /* Any eastward offset to add while reading */
  float north;  /* Any northward offset to add while reading */
  int docont=0; /* Default to read into the normal model */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("rmodel", OB_INDEX))
    return -1;
/*
 * Get arguments.
 */
  switch(npar) {
  case 2:
    docont = *LOGPTR(invals[1]);
  case 1: default:
    modfil = *STRPTR(invals[0]);
  };
/*
 * Check that the requested model file exists before initializing.
 */
  if(!file_exists(modfil)) {
    lprintf(stderr, "rmodel: File \"%s\" does not exist\n", modfil);
    return -1;
  };
/*
 * The residual map will have to be re-made before use.
 */
  if(vlbmap)
    vlbmap->domap = MAP_IS_STALE;
/*
 * Clear any existing models of the same type.
 */
  if(docont)
    clrmod(vlbob, 0, 0, 1);
  else
    clrmod(vlbob, 1, 1, 0);
/*
 * Determine the coordinate offsets to add to the model components.
 */
  east  = vlbob->geom.east * rtomas;
  north = vlbob->geom.north * rtomas;
/*
 * Read the new model.
 */
  if(rmodel(docont ? vlbob->cnewmod:vlbob->newmod, east, north, 1, modfil))
    return -1;
  return no_error;
}

/*.......................................................................
 * Apply amplitude self-cal to find and apply overall amplitude factors for
 * each telescope over the whole observation.
 *
 * Input:
 *   dofloat  char  If TRUE then allow the amplitudes to float without
 *                  restraint. Otherwise normalise amplitude corrections.
 */
static Template(gscal_fn)
{
  int dofloat=0;  /* Prevent unconstrained amplitude corrections by default */
  int doamp=1;    /* Default to doing amplitude corrections */
  int dophs=0;    /* Default to not doing phase corrections */
  int flagged;    /* Returned true by slfcal if any data were flagged. */
  int iret;       /* Return value from slfcal */
  if(npar > 0)
    dofloat = *LOGPTR(invals[0]);
/*
 * Make sure that an observation has been read.
 */
  if(nodata("gscale", OB_SELECT))
    return -1;
/*
 * Mark the map as invalid since self-cal changes the visibilities.
 * Also mark the beam as invalid if error weighting is currently selected,
 * because amplitude corrections are also applied to the amplitude
 * uncertainties.
 */
  if(vlbmap) {
    vlbmap->domap = MAP_IS_STALE;
    if(invpar.errpow < 0.0)
      vlbmap->dobeam = 1;
  };
/*
 * Amplitude self-cal for overall amplitude errors.
 */
  lprintf(stdout, "Performing overall amplitude self-cal\n");
  iret = slfcal(vlbob, -1, 1, slfpar.gauval, slfpar.gaurad,
		0.0, doamp, dophs, dofloat,
		(doamp ? slfpar.a_mintel : slfpar.p_mintel),
		slfpar.doflag, 1, slfpar.maxamp, slfpar.maxphs,
		invpar.uvmin, invpar.uvmax, &flagged);
/*
 * If data were flagged then mark the beam as invalid.
 */
  if(vlbmap && flagged)
    vlbmap->dobeam=1;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Establish the latest CLEAN model by appending it to the current
 * established model and adding it to the model visiblities in the UV
 * plane.
 */
static Template(keep_fn)
{
/*
 * Make sure that an observation has been read.
 */
  if(nodata("keep", OB_SELECT))
    return -1;
/*
 * Make sure that there are components to be established, so that the
 * map is not marked as invalid unnecessarily.
 */
  if(vlbob->newmod->ncmp + vlbob->cnewmod->ncmp > 0) {
/*
 * The UV data will have to be re-inverted to produce the modified
 * residual map. (No need to flag the beam for re-calculation).
 */
    if(vlbmap)
      vlbmap->domap = MAP_IS_STALE;
/*
 * Establish the tentative model.
 */
    if(mergemod(vlbob, 1))
      return -1;
  };
  return no_error;
}

/*.......................................................................
 * Clear both the model visibilities and the two model component lists.
 *
 * Input:
 *  doold   bool  If TRUE clear the established model.
 *  donew   bool  If TRUE clear the tentative model.
 *  docont  bool  If TRUE clear the continuum model.
 */
static Template(clrmod_fn)
{
  int doold=0;   /* Default to not clearing the established model */
  int donew=1;   /* Default to clearing the tentative model */
  int docont=0;  /* Default to not clearing the continuum model */
/*
 * Check that we have an observation.
 */
  if(nodata("clrmod", OB_INDEX))
    return 1;
/*
 * Override defaults if user gave an argument.
 */
  switch(npar) {
  case 3:
    docont = *LOGPTR(invals[2]);
  case 2:
    donew = *LOGPTR(invals[1]);
  case 1:
    doold = *LOGPTR(invals[0]);
  };
/*
 * Clear the model(s).
 */
  clrmod(vlbob, doold, donew, docont);
/*
 * The modified residual map will have to be recalculated before its
 * next use. (The beam is unchanged).
 */
  if(vlbmap)
    vlbmap->domap = MAP_IS_STALE;
  return no_error;
}

/*.......................................................................
 * Modify the UV data and established model such that the map is shifted
 * by xshift eastwards and yshift Northwards. Also shift any clean
 * windows that are currently defined.
 *
 * Input:
 *  xshift  float  The distance to shift Eastwards (map xy units).
 *  yshift  float  The distance to shift Northwards (map xy units).
 */
static Template(shift_fn)
{
/*
 * Check that there is something to shift.
 */
  if(nodata("shift", OB_INDEX))
    return -1;
/*
 * Are we changing the current shift?
 */
  if(npar > 0) {
    float xshift = 0.0f;
    float yshift = 0.0f;
/*
 * Get the required shifts (converted to radians).
 */
    switch(npar) {
    default:
    case 2:
      yshift = xytorad(*FLTPTR(invals[1]));
    case 1:
      xshift = xytorad(*FLTPTR(invals[0]));
    };
/*
 * Report the shifts before applying them.
 */
    lprintf(stdout,
	    "Shifting UV data, models and windows by: %g (%s) East\n",
	    radtoxy(xshift), mapunits(U_NAME));
    lprintf(stdout,
	    "Shifting UV data, models and windows by: %g (%s) North\n",
	    radtoxy(yshift), mapunits(U_NAME));
/*
 * First shift the clean windows.
 */
    shiftwin(vlbwins, xshift, yshift);
/*
 * Shift the observed visibilities and models.
 */
    if(obshift(vlbob, xshift, yshift))
      return -1;
/*
 * The shifted map will need to be recalculated before its next use.
 * (The beam is unchanged).
 */
    if(vlbmap)
      vlbmap->domap = MAP_IS_STALE;
  };
/*
 * Report the accumulated total shifts.
 */
  lprintf(stdout, "Total accumulated eastward shift  = %g (%s).\n",
	  radtoxy(vlbob->geom.east), mapunits(U_NAME));
  lprintf(stdout, "Total accumulated northward shift = %g (%s).\n",
	  radtoxy(vlbob->geom.north), mapunits(U_NAME));
  return no_error;
}

/*.......................................................................
 * Remove any accumulated position shifts applied to the phase center,
 * CLEAN windows and models.
 */
static Template(unshift_fn)
{
/*
 * Check that there something to unshift.
 */
  if(nodata("unshift", OB_INDEX))
    return -1;
/*
 * Report action.
 */
  lprintf(stdout, "unshift: Removing accumulated position shifts.\n");
/*
 * First shift the clean windows back.
 */
  shiftwin(vlbwins, -vlbob->geom.east, -vlbob->geom.north);
/*
 * Shift the observed visibilities and established model.
 */
  if(obunshift(vlbob))
    return -1;
/*
 * The shifted map will need to be recalculated before its next use.
 * (The beam is unchanged).
 */
  if(vlbmap)
    vlbmap->domap = MAP_IS_STALE;
  return no_error;
}

/*.......................................................................
 * Average an observation and remove the resulting unused integrations.
 *
 * Input:
 *  av_time  float  The averaging time (seconds).
 *  doscat    char  If true estimate amplitude/phase uncertainties from
 *                  the scatter of the input data. (default=no scatter).
 */
static Template(uvav_fn)
{
  float av_time=0.0f; /* The required averaging time (seconds) */
  int doscat=0;       /* No scatter estimates by default */
/*
 * Check that we have a UV data-set to be averaged.
 */
  if(nodata("uvaver", OB_INDEX))
    return -1;
/*
 * Get optional arguments. Case fallthrough is deliberate.
 */
  switch(npar) {
  case 2:
    doscat = *LOGPTR(invals[1]);
  case 1:
    av_time = *FLTPTR(invals[0]);
  };
  if(av_time <= 0.0f) {
    lprintf(stderr, "uvaver: Illegal averaging time (%g)\n",av_time);
    return -1;
  };
/*
 * Average.
 */
  vlbob=uvaver(vlbob, av_time, doscat);
/*
 * Both the map and beam wil have to be calculated before their next use.
 */
  if(vlbmap)
    vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
  return no_error;
}

/*.......................................................................
 * List useful parameters from an Observation header.
 */
static Template(head_fn)
{
/*
 * Check that we have a UV data-set to describe.
 */
  if(nodata("header", OB_INDEX))
    return -1;
/*
 * List the header parameters.
 */
  vlbhead(vlbob);
  return no_error;
}

/*.......................................................................
 * Undo recorded telescope corrections.
 *
 * Input:
 *  dophs   logical  Undo telescope phase corrections if true.
 *                   Default=false.
 *  doamp   logical  Undo telescope amplitude corrections if true.
 *                   Default=false.
 *  doflag  logical  Undo telescope correction flags if true.
 *                   Default=false.
 */
static Template(uncal_fn)
{
  int dophs=0;  /* Default to not undoing of phase corrections */
  int doamp=0;  /* Default to not undoing of amplitude corrections */
  int doflag=0; /* Default to not undoing of correction flags */
/*
 * Check that we have a UV data-set to uncorrect.
 */
  if(nodata("uncalib", OB_INDEX))
    return -1;
/*
 * Override defaults with user arguments (Note fallthrough).
 */
  switch(npar) {
  case 3:
    doflag = *LOGPTR(invals[2]);
  case 2:
    doamp = *LOGPTR(invals[1]);
  case 1:
    dophs = *LOGPTR(invals[0]);
  };
/*
 * NULL operation?
 */
  if(!doamp && !dophs && !doflag) {
    return no_error;
  } else {
    if(vlbmap)
      vlbmap->domap=MAP_IS_STALE;
/*
 * Undo telescope corrections.
 */
    uncalib(vlbob, doamp, dophs, doflag, 1);
/*
 * Keep user informed about what has been done.
 */
    if(dophs)
      lprintf(stdout,
	"uncal: All telescope phase corrections have been un-done.\n");
    if(doamp)
      lprintf(stdout,
        "uncal: All telescope amplitude corrections have been un-done.\n");
    if(doflag)
      lprintf(stdout,
        "uncal: All telescope correction flags have been un-done.\n");
  };
  return no_error;
}

/*.......................................................................
 * Plot the amplitude and phase corrections for a given station.
 *
 * Input:
 *  tname  char *  Name of telescope who's corrections should be plotted.
 *                 (Optional) If omitted interactive cursor mode will be
 *                 invoked.
 *  cif     int    The start IF, or 0 for the default.
 */
static Template(corpl_fn)
{
  Telspec *ts=NULL;/* Telescope specification string */
  int cif = -1;    /* Default start IF index */
  int docurs=1;    /* Request interactive cursor mode if possible */
  int modified=0;  /* Returned true by corplot() if the data were modified */
  int iret;        /* Return value of corplot() */
/*
 * Check that we have a UV data-set to plot.
 */
  if(nodata("corplot", OB_INDEX))
    return -1;
/*
 * Get the telescope number and sub-array to be plotted.
 */
  switch(npar) {
  case 2:
    cif = *INTPTR(invals[1]) - 1;
  case 1:
    ts = read_Telspec(vlbob, *STRPTR(invals[0]), NULL, 0);
    if(!ts)
      return -1;
  };
/*
 * Ensure that PGPLOT is open.
 */
  if(make_open() == -1)
    return -1;
/*
 * Display the corrections.
 */
  iret = corplot(vlbob, ts, cif, docurs, &modified);
/*
 * Edits to the corrections will invalidate the map and beam.
 */
  if(vlbmap && modified)
    vlbmap->domap = vlbmap->dobeam = MAP_IS_STALE;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Return the name of a telescope given its index.
 *
 * Input:
 *   itel    int    The index of the telescope.
 *   isub    int    The index of the sub-array - defaults to 1.
 * Output:
 *   tname  char *  The name of telescope 'itel'.
 */
static Template(tname_fn)
{
  Subarray *sub; /* The descriptor of the requested sub-array */
  int itel=0;    /* Requested telescope index */
  int isub=0;    /* Requested sub-array index */
  size_t slen;   /* The length of the telescope name string */
  char *cptr;
/*
 * Check that we have a UV data-set to examine.
 */
  if(nodata("telname", OB_INDEX))
    return -1;
/*
 * Get the user arguments.
 */
  switch(npar) {
  case 2:
    isub = *INTPTR(invals[1]) - 1;
  case 1:
    itel = *INTPTR(invals[0]) - 1;
  };
/*
 * Valid sub-array index?
 */
  if(isub<0 || isub>=vlbob->nsub) {
    lprintf(stderr, "telname: Out of range sub-array index: %d\n", isub+1);
    return -1;
  };
/*
 * Get the descriptor of the requested sub-array.
 */
  sub = &vlbob->sub[isub];
/*
 * Check the telescope index.
 */
  if(itel < 0 || itel >= sub->nstat) {
    lprintf(stderr, "telname: Out of range telescope index: %d\n", itel);
    return -1;
  };
/*
 * Determine the length of the required return string, sufficient to
 * hold the sub-array number (assume up to 3 digits) and the telescope
 * name.
 */
  slen = strlen(sub->tel[itel].name) + 4;
/*
 * Now allocate memory for a copy of the telescope name string.
 */
  cptr = stralloc(slen * sizeof(char));
  if(cptr==NULL) {
    lprintf(stderr, "telname: Insufficient memory for copy of name.\n");
    return -1;
  };
/*
 * Copy the telescope name into the new string and assign it to its return
 * slot.
 */
  sprintf(cptr, "%d:%s", isub+1, sub->tel[itel].name);
  *STRPTR(outvals) = cptr;
  return no_error;
}

/*.......................................................................
 * Return the count of the number of telescopes in the observation.
 *
 * Input:
 *   isub    int    The index of the sub-array - defaults to 1.
 * Output:
 *   return  int    The number of telescopes in sub-array isub.
 */
static Template(ntel_fn)
{
  Subarray *sub; /* The descriptor of the requested sub-array */
  int isub=0;    /* The index of the sub-array */
/*
 * Check that we have a UV data-set to examine.
 */
  if(nodata("ntel", OB_INDEX))
    return -1;
/*
 * Get the user arguments.
 */
  switch(npar) {
  case 1:
    isub = *INTPTR(invals[0]) - 1;
  };
/*
 * Valid sub-array index?
 */
  if(isub<0 || isub>=vlbob->nsub) {
    lprintf(stderr, "ntel: Out of range sub-array index: %d\n", isub+1);
    return -1;
  };
/*
 * Get the descriptor of the requested sub-array.
 */
  sub = &vlbob->sub[isub];
/*
 * Assign the number of telescopes in the sub-array as the return
 * value.
 */
  *INTPTR(outvals) = sub->nstat;
  return no_error;
}

/*.......................................................................
 * Return the name of a baseline given its index.
 *
 * Input:
 *   ibase   int    The index of the baseline.
 *   isub    int    The index of the sub-array - defaults to 1.
 * Output:
 *   bname  char *  The name of baseline 'ibase'.
 */
static Template(bname_fn)
{
  Subarray *sub;/* The descriptor of the requested sub-array */
  int ibase=0;  /* The index of the requested baseline */
  int isub=0;   /* The index of the requested sub-array */
  char *name1;  /* Name of first telescope on baseline */
  char *name2;  /* Name of second telescope on baseline */
  int slen;     /* Length of baseline string */
  char *cptr;
/*
 * Check that we have a UV data-set to examine.
 */
  if(nodata("basename", OB_INDEX))
    return -1;
/*
 * Get the user arguments.
 */
  switch(npar) {
  case 2:
    isub = *INTPTR(invals[1]) - 1;
  case 1:
    ibase = *INTPTR(invals[0]) - 1;
  };
/*
 * Valid sub-array index?
 */
  if(isub<0 || isub>=vlbob->nsub) {
    lprintf(stderr, "basename: Out of range sub-array index: %d\n", isub+1);
    return -1;
  };
/*
 * Get the descriptor of the requested sub-array.
 */
  sub = &vlbob->sub[isub];
/*
 * Check the baseline index.
 */
  if(ibase < 0 || ibase >= sub->nbase) {
    lprintf(stderr, "basename: Out of range baseline index: %d\n", ibase);
    return -1;
  };
/*
 * Get pointers to the two telescope name strings.
 */
  name1 = sub->tel[sub->base[ibase].tel_a].name;
  name2 = sub->tel[sub->base[ibase].tel_b].name;
/*
 * See how long a string will be required to hold a string combining
 * the sub-array number (assume up to 3 digits) and the two telescope
 * names separated by a hyphen. (Don't worry about the
 * '\0' on the end of the combined string - stralloc handles that).
 */
  slen = strlen(name1) + strlen(name2) + 4;
/*
 * Now allocate memory for the baseline string.
 */
  cptr = stralloc(slen);
  if(cptr==NULL) {
    lprintf(stderr, "basename: Insufficient memory for baseline string\n");
    return -1;
  };
/*
 * Compile the baseline string.
 */
  sprintf(cptr,"%d:%s-%s", isub+1, name1, name2);
/*
 * Copy it to its output slot.
 */
  *STRPTR(outvals) = cptr;
  return no_error;
}

/*.......................................................................
 * Return the count of the number of baselines in the observation.
 *
 * Input:
 *  isub    int    The index of the sub-array - defaults to 1.
 * Output:
 *  return  int    The number of baselines in the sub-array.
 */
static Template(nbase_fn)
{
  Subarray *sub; /* The descriptor of the requested sub-array */
  int isub=0;    /* The index of the sub-array */
/*
 * Check that we have a UV data-set to plot.
 */
  if(nodata("nbase", OB_INDEX))
    return -1;
/*
 * Get the user arguments.
 */
  switch(npar) {
  case 1:
    isub = *INTPTR(invals[0]) - 1;
  };
/*
 * Valid sub-array index?
 */
  if(isub<0 || isub>=vlbob->nsub) {
    lprintf(stderr, "nbase: Out of range sub-array index: %d\n", isub+1);
    return -1;
  };
/*
 * Get the descriptor of the requested sub-array.
 */
  sub = &vlbob->sub[isub];
/*
 * Assign the number of telescopes in the sub-array as the return
 * value.
 */
  *INTPTR(outvals) = sub->nbase;
  return no_error;
}

/*.......................................................................
 * Return the count of the number of sub-arrays in the observation.
 *
 * Output:
 *  return  int    The number of sub-arrays in the observation.
 */
static Template(nsub_fn)
{
/*
 * Check that we have a UV data-set to plot.
 */
  if(nodata("nsub", OB_INDEX))
    return -1;
/*
 * Assign the number of sub-arrays as the return value.
 */
  *INTPTR(outvals) = vlbob->nsub;
  return no_error;
}

/*.......................................................................
 * Return the count of the number of IFs in the observation.
 *
 * Output:
 *  return  int    The number of IFs in the observation.
 */
static Template(nif_fn)
{
/*
 * Check that we have a UV data-set to plot.
 */
  if(nodata("nif", OB_INDEX))
    return -1;
/*
 * Assign the number of IFs as the return value.
 */
  *INTPTR(outvals) = vlbob->nif;
  return no_error;
}

/*.......................................................................
 * Return the count of the number of channels in the observation.
 *
 * Output:
 *  return  int    The number of spectral-line channels in the observation.
 */
static Template(nchan_fn)
{
/*
 * Check that we have a UV data-set to plot.
 */
  if(nodata("nchan", OB_INDEX))
    return -1;
/*
 * Assign the number of channels as the return value.
 */
  *INTPTR(outvals) = vlbob->nchan;
  return no_error;
}

/*.......................................................................
 * Add a new clean window to the list of windows.
 *
 * Input:
 *  xa  float  One edge of the window along the X axis.
 *  xb  float  The other edge of the window along the X axis.
 *  ya  float  One edge of the window along the Y axis.
 *  yb  float  The other edge of the window along the Y axis.
 */
static Template(addwin_fn)
{
  float xa,xb,ya,yb;
/*
 * Either no arguments or 4 arguments are required.
 */
  if(npar != 4) {
    lprintf(stderr, "addwin: Insufficient arguments\n");
    return -1;
  };
/*
 * Get the arguments.
 */
  xa = xytorad(*FLTPTR(invals[0]));
  xb = xytorad(*FLTPTR(invals[1]));
  ya = xytorad(*FLTPTR(invals[2]));
  yb = xytorad(*FLTPTR(invals[3]));
/*
 * First create an empty Mapwin list if one hasn't been created yet.
 */
  if(vlbwins==NULL && (vlbwins=new_Mapwin()) == NULL)
    return -1;
/*
 * Add the new window.
 */
  if(add_win(vlbwins, xa, xb, ya, yb)==NULL)
    return -1;
  return no_error;
}

/*.......................................................................
 * Delete all clean windows.
 */
static Template(delwin_fn)
{
  vlbwins=del_Mapwin(vlbwins);
  lprintf(stdout, "All clean windows deleted\n");
  return no_error;
}

/*.......................................................................
 * Delete all components that lie inside or (optionally) outside the
 * current clean windows.
 *
 * Input:
 *  doout  logical  (optional - default=false). If true then delete
 *                  all components that lie inside the current clean
 *                  windows. If false then delete those outside all
 *                  windows.
 */
static Template(winmod_fn)
{
  int doout;      /* If true, keep the components outside the windows */
/*
 * No data?
 */
  if(nodata("winmod", OB_INDEX))
    return -1;
/*
 * Interpret user arguments.
 */
  doout = (npar>0) ? *LOGPTR(invals[0]) : 0;
/*
 * No clean windows to apply?
 */
  if(vlbwins==NULL) {
    lprintf(stderr, "winmod: There are no clean windows\n");
    return -1;
  };
/*
 * The map will have to be recomputed after modifying the model.
 */
  if(vlbmap != NULL)
    vlbmap->domap=MAP_IS_STALE;
/*
 * Window the established and tentative models.
 */
  if(obwinmod(vlbob, vlbwins, doout))
    return -1;
  return no_error;
}

/*.......................................................................
 * Read a starting model from a given model file, perform phase self-cal
 * using it and then delete the model.
 *
 * Input:
 *   modfile  char *  The name of the model file.
 */
static Template(startmod_fn)
{
/*
 * Declare descriptor arguments for self_fn{}.
 */
  static char doamp = '\0';
  static char dofloat = '\0';
  static float solint = 0.0f;
  static Descriptor d_doamp  = {'l',0,R_ONLY,1,{1,1,1},&doamp};
  static Descriptor d_dofloat = {'l',0,R_ONLY,1,{1,1,1},&dofloat};
  static Descriptor d_solint  = {'f',0,R_ONLY,1,{1,1,1},&solint};
  static Descriptor *self_args[3] = {&d_doamp, &d_dofloat, &d_solint};
/*
 * Is there anything to be self-cal'd?
 */
  if(nodata("startmod", OB_SELECT))
    return -1;
/*
 * Clear the existing models.
 */
  clrmod(vlbob, 1, 1, 1);
/*
 * If a model file name was given, open the file and read the starting model.
 */
  if(npar>0 && *STRPTR(invals[0])[0] != '\0') {
    if(rmodel_fn(invals,1,outvals)==-1)
      return -1;
  }
/*
 * If no model file was named, create a single component model.
 */
  else {
    lprintf(stdout, "Applying default point source starting model.\n");
    if(add_xycmp(vlbob->newmod, 1,0, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, M_DELT,
		 0.0, 0.0) == NULL) {
      return -1;
    };
  };
/*
 * Set up self_fn arguments.
 */
  *LOGPTR(&d_doamp) = 0;
  *LOGPTR(&d_dofloat) = 0;
  *FLTPTR(&d_solint) = npar > 1 ? *FLTPTR(invals[1]) : 0.0f;
/*
 * Use the new model to perform phase self-cal.
 */
  if(self_fn(self_args, 3, outvals) == -1)
    return -1;
/* 
 * We have now finished with the model so delete it.
 */
  clrmod(vlbob, 1, 1, 1);
  lprintf(stdout, "Redundant starting model cleared.\n");
  return no_error;
}

/*.......................................................................
 * Plot amplitude versus UV radius.
 *
 * Input:
 *  telnam  char * The name of a telescope to highlight.
 *  uvmin  float   The min UV radius to display (user UV distance units).
 *  uvmax  float   The max UV radius (or 0 for full range) to plot.
 *                 (User UV distance units).
 *  ampmin float   The min amplitude to be displayed.
 *  ampmax float   The max amplitude to be displayed.
 *  phsmin float   The min phase to be displayed (degrees).
 *  phsmax float   The max phase to be displayed (degrees).
 *  docur logical  Default TRUE. If false, bypass interactive mode.
 */
static Template(uvrad_fn)
{
  Telspec *ts=NULL;  /* Telescope specification descriptor */
  char *opts="m1";   /* Default plot options */
  int docur=1;       /* Default to use of cursor where available */
  float uvmin=0.0f;  /* Default to full UV range */
  float uvmax=0.0f;
  float ampmin=0.0f; /* Default to full amplitude range */
  float ampmax=0.0f;
  float phsmin=0.0f; /* Default to full phase range */
  float phsmax=0.0f;
  int modified;      /* Returned true by uvradplt() if data were modified */
  int iret;          /* Return value from uvradplt() */
/*
 * Make sure that we have data to be plotted.
 */
  if(nodata("radplot", OB_SELECT))
    return -1;
/*
 * Selectively overide defaults where given. NB. Fallthrough is deliberate.
 */
  switch(npar) {
  case 8:
    docur = *LOGPTR(invals[7]);
  case 7:
    phsmax = *FLTPTR(invals[6]) * dtor;
  case 6:
    phsmin = *FLTPTR(invals[5]) * dtor;
  case 5:
    ampmax = *FLTPTR(invals[4]);
  case 4:
    ampmin = *FLTPTR(invals[3]);
  case 3:
    uvmax = uvtowav(*FLTPTR(invals[2]));
  case 2:
    uvmin = uvtowav(*FLTPTR(invals[1]));
  case 1:
    ts = read_Telspec(vlbob, *STRPTR(invals[0]), NULL, 0);
    if(!ts)
      return -1;
  };
/*
 * Get the start telescope and sub-array indexes.
 */
/*
 * Establish the latest clean model in order to be able to display it.
 */
  if(keep_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Open the display device.
 */
  if(make_open() == -1)
    return -1;
/*
 * Does the rflags options string contain any entries?
 * If so override the default options string.
 */
  if(**STRPTR(&rflags)=='\0')
    lprintf(stdout, "Using default options string \"%s\"\n", opts);
  else {
    opts = *STRPTR(&rflags);
    lprintf(stdout,
       "Overriding default options with user defined rflags=\"%s\"\n", opts);
  };
/*
 * Plot the data.
 */
  iret = uvradplt(vlbob, ts, docur, opts, 0, 0.0f, uvmin, uvmax,
		  ampmin, ampmax, phsmin, phsmax, &modified);
/*
 * If the data were modified, mark the map and beam as out of date.
 */
  if(vlbmap && modified)
    vlbmap->domap = vlbmap->dobeam = MAP_IS_STALE;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Apply phase (and amplitude) self-cal to the UV data in vlbob wrt the
 * established and tentative models.
 *
 * Input:
 *   doamp    char  (Default=FALSE) If TRUE then do solve for amplitude
 *                  corrections in addition to phases.
 *   dofloat  char  (Default=FALSE). If TRUE then allow the amplitudes
 *                  to float without restraint. Otherwise normalise
 *                  amplitude corrections.
 *   solint  float  The solution interval in minutes.
 */
static Template(self_fn)
{
  int doamp=0;    /* Prevent amplitude corrections by default */
  int dophs=1;    /* Perform phase self-cal by default */
  int dofloat=0;  /* Prevent unconstrained amplitude corrections by default */
  float solint=0.0f;/* Default to integration-by-integration solutions */
  int flagged;    /* Returned true by slfcal if any data were flagged */
  int iret;       /* Return value from slfcal() */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("selfcal", OB_SELECT))
    return -1;
/*
 * Overide defualts with arguments where given. Note the intentional
 * fallthrough in the switch.
 */
  switch (npar) {
  case 3:
    solint = *FLTPTR(invals[2]);
  case 2:
    dofloat = *LOGPTR(invals[1]);
  case 1:
    doamp = *LOGPTR(invals[0]);
  };
/*
 * Mark the map as invalid since self-cal changes the visibilities.
 * If amplitude self-cal has been requested, mark the beam as invalid
 * if error weighting is currently selected since amplitude corrections
 * are also applied to the amplitude uncertainties.
 */
  if(vlbmap) {
    vlbmap->domap = MAP_IS_STALE;
    if(doamp && invpar.errpow < 0.0)
      vlbmap->dobeam = 1;
  };
/*
 * Keep the user informed.
 */
  lprintf(stdout, "Performing %s self-cal", doamp?"amp+phase":"phase");
  if(solint > 0.0f)
    lprintf(stdout, " over %g minute time intervals\n", solint);
  else
    lprintf(stdout, "\n");
/*
 * Apply self-cal.
 */
  iret = slfcal(vlbob, -1, 1, slfpar.gauval, slfpar.gaurad,
		solint, doamp, dophs, dofloat,
		(doamp ? slfpar.a_mintel : slfpar.p_mintel),
		slfpar.doflag, 0, slfpar.maxamp, slfpar.maxphs,
		invpar.uvmin, invpar.uvmax, &flagged);
/*
 * If data were flagged then mark the beam as invalid.
 */
  if(vlbmap && flagged)
    vlbmap->dobeam=1;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Specify a UV gaussian taper for weighting down short baselines in
 * self-cal.
 * If the specified gaussian value is not between 0 and 1 or the
 * specified radius is <= 0 then cancel any previous taper.
 *
 * Input:
 *  gauval    float  The gaussian taper value at 'gaurad'.
 *  gaurad    float  The radius (user UV distance units) that the taper has
 *                   value 'gauval'.
 */
static Template(staper_fn)
{
/*
 * The taper is undefined until an observation has been read.
 */
  if(nodata("selftaper", OB_INDEX))
    return -1;
/*
 * Are we changing the existing taper state?
 */
  if(npar > 0) {
/*
 * Override current taper values.
 */
    switch(npar) {
    case 2:
      slfpar.gaurad = uvtowav(*FLTPTR(invals[1]));  /* Convert to wavelengths */
    case 1:
      slfpar.gauval = *FLTPTR(invals[0]);
    };
/*
 * Cancel the taper if the new values are out of bounds.
 */
    if(slfpar.gauval<=0.0f || slfpar.gauval>=0.99f || slfpar.gaurad<=0.0f) {
      slfpar.gauval = 0.0f;
      slfpar.gaurad = 0.0f;
    };
  };
/*
 * Report the current settings.
 */
  if(slfpar.gauval <= 0.0f || slfpar.gaurad <= 0.0f) {
    lprintf(stdout, "No selfcal UV-taper is currently set.\n");
  } else {
    lprintf(stdout, "Selfcal gaussian taper: value %g at UV radius = %g %s.\n",
	    slfpar.gauval, wavtouv(slfpar.gaurad), uvwunits(U_TLAB));
  };
  return no_error;
}

/*.......................................................................
 * Assign or cancel limits to self-cal amplitude and phase corrections.
 *
 * Input:
 *  maxamp  float  The max amplitude correction allowed (>1.0).
 *  maxphs  float  The max phase correction allowed (>0.0 degrees).
 */
static Template(slims_fn)
{
/*
 * The taper is undefined until an observation has been read.
 */
  if(nodata("selflims", OB_INDEX))
    return -1;
/*
 * Override current values where given.
 */
  if(npar > 0) {
    switch (npar) {
    default:
    case 2:
      slfpar.maxphs = *FLTPTR(invals[1]) * dtor; /* Phase in radians */
    case 1:
      slfpar.maxamp = *FLTPTR(invals[0]);
    };
/*
 * Enforce limits.
 */
    if(slfpar.maxphs < 0.0f || slfpar.maxphs >= 180.0f)
      slfpar.maxphs = 0.0f;
    if(slfpar.maxamp <= 1.0f)
      slfpar.maxamp = 0.0f;
  };
/*
 * Report current amplitude limits.
 */
  if(slfpar.maxamp > 0.0f)
    lprintf(stdout, "Selfcal amplitude corrections will be limited to %g -> %g.\n", 1.0f/slfpar.maxamp,
	    slfpar.maxamp);
  else
    lprintf(stdout, "Selfcal amplitude corrections will not be limited.\n");
/*
 * Report new phase limits.
 */
  if(slfpar.maxphs > 0.0f)
    lprintf(stdout, "Selfcal phase corrections will be limited to %g -> %g degrees.\n",
	   -slfpar.maxphs * rtod, slfpar.maxphs * rtod);
  else
    lprintf(stdout, "Selfcal phase corrections will not be limited.\n");
  return no_error;
}

/*.......................................................................
 * Read a new clean-window list from a given file. Any existing windows
 * will be deleted.
 *
 * Input:
 *  winfile  char *  The name of the file containing the windows.
 */
static Template(rwins_fn)
{
  char *winfile = *STRPTR(invals[0]);
  float xshift;  /* Any eastward offset to add while reading */
  float yshift;  /* Any northward offset to add while reading */
/*
 * Check that the requested window file exists before initializing.
 */
  if(!file_exists(winfile)) {
    lprintf(stderr, "rwins: File \"%s\" does not exist\n", winfile);
    return -1;
  };
/*
 * Delete any existing windows.
 */
  vlbwins = del_Mapwin(vlbwins);
/*
 * Determine the coordinate offsets to add to the windows.
 */
  xshift = vlbob ? vlbob->geom.east : 0.0f;
  yshift = vlbob ? vlbob->geom.north : 0.0f;
/*
 * Create a new empty list of windows.
 */
  vlbwins = new_Mapwin();
/*
 * Read the windows from the given file.
 */
  if(vlbwins==NULL || rwins(vlbwins, winfile, xshift, yshift))
    return -1;
  return no_error;
}

/*.......................................................................
 * Write the current clean windows to a file.
 *
 * Input:
 *  winfile  char *  The name of the file to contain the windows.
 */
static Template(wwins_fn)
{
  int do_old=0; /* By default use new rationalized window format */
  char *winfile=NULL;
  float xshift;  /* Any eastward offset to remove while writing */
  float yshift;  /* Any northward offset to remove while writing */
/*
 * Any windows to write?
 */
  if(vlbwins==NULL || vlbwins->nwin==0) {
    lprintf(stderr, "wwins: No CLEAN windows to write.\n");
    return -1;
  };
/*
 * Get user arguments.
 */
  switch(npar) {
  case 2:
    do_old = *LOGPTR(invals[1]);
  case 1:
    winfile = *STRPTR(invals[0]);
  };
/*
 * Determine the coordinate offsets to remove while writing.
 */
  xshift = vlbob ? vlbob->geom.east : 0.0f;
  yshift = vlbob ? vlbob->geom.north : 0.0f;
/*
 * Write the windows.
 */
  if(wwins(vlbwins, winfile, xshift, yshift, do_old) != 0)
    return -1;
  return no_error;
}

/*.......................................................................
 * Plot visibilities with the option of cursor control.
 *
 * Input:
 *  nrow     int    The initial number of rows in the plot.
 *                  Default=0 which gives ob->nstat-1.
 *  basestr char *  The initial baseline specification.
 *  cif      int    The IF to plot first.
 *  npage    int    The max number of pages to plot when non-interactive.
 */
static Template(vplot_fn)
{
  Basespec *bs=NULL; /* Baseline specification descriptor */
  int nrow=0;        /* Default number of rows */
  int modified=0;    /* Returned true by vedit if the data were modified */
  int cif = -1;      /* The start IF index */
  int iret;          /* Return code from vedit() */
  int npage = 0;     /* The max number of pages to plot when non-interactive */
  char *opts="efbm3";/* Default plot options */
/*
 * Make sure that we have data to be plotted.
 */
  if(nodata("vplot", OB_SELECT))
    return -1;
/*
 * Make sure that a PGPLOT device is open.
 */
  if(make_open() == -1)
    return -1;
/*
 * Check arguments.
 */
  switch (npar) {
  case 4:
    npage = *INTPTR(invals[3]);
  case 3:
    cif = *INTPTR(invals[2]) - 1;
  case 2:
    bs = read_Basespec(vlbob, *STRPTR(invals[1]), NULL, 0);
    if(!bs)
      return -1;
  case 1:
    nrow = *INTPTR(invals[0]);
  };
/*
 * Does the vflags options string contain any entries?
 * If so override the default options string.
 */
  if(**STRPTR(&vflags)=='\0')
    lprintf(stdout, "Using default options string \"%s\"\n", opts);
  else {
    opts = *STRPTR(&vflags);
    lprintf(stdout,
       "Overriding default options with user defined vflags=\"%s\"\n", opts);
  };
/*
 * Establish the latest clean model so that we have an up-to-date
 * representation of the UV model to display.
 */
  if(keep_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Plot the data.
 */
  iret = vedit(vlbob, bs, cif, nrow, npage, 1, opts, 0, 1, 1, 0, 0, 0, 0,
	       &modified);
/*
 * Newly flagged or unflagged visibilities affect the map and beam.
 * If vedit() reported the data as modified, mark both the map and beam
 * as out of date.
 */
  if(vlbmap && modified)
    vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Plot U versus V-distance.
 *
 * Input:
 *  telnam  char * The name of a telescope to highlight.
 *  umax   float   The min UV radius to display (User UV distance units).
 *  vmax   float   The max UV radius (or 0 for full range) to plot.
 *                 (user UV distance units).
 *  docur logical  Default TRUE. If false, bypass interactive mode.
 */
static Template(uvplt_fn)
{
  Telspec *ts=NULL;  /* Telescope specification descriptor */
  int docur=1;       /* Default to use of cursor where available */
  float umax=0.0f;   /* Default to full UV range */
  float vmax=0.0f;
  int modified;      /* Returned true by uvplot() if data were modified */
  int iret;          /* Return value from uvplot() */
  char *opts="";     /* Default plot options */
/*
 * Make sure that we have data to be plotted.
 */
  if(nodata("uvplot", OB_SELECT))
    return -1;
/*
 * Selectively overide defaults where given. NB. Fallthrough is deliberate.
 */
  switch(npar) {
  case 4:
    docur = *LOGPTR(invals[3]);
  case 3:
    vmax = uvtowav(*FLTPTR(invals[2]));
  case 2:
    umax = uvtowav(*FLTPTR(invals[1]));
  case 1:
    ts = read_Telspec(vlbob, *STRPTR(invals[0]), NULL, 0);
    if(!ts)
      return -1;
  };
/*
 * Open the display device.
 */
  if(make_open() == -1)
    return -1;
/*
 * Does the tflags options string contain any entries?
 * If so override the default options string.
 */
  if(**STRPTR(&uflags)=='\0')
    lprintf(stdout, "Using default options string \"%s\"\n", opts);
  else {
    opts = *STRPTR(&uflags);
    lprintf(stdout,
       "Overriding default options with user defined uflags=\"%s\"\n", opts);
  };
/*
 * Plot the data.
 */
  iret = uvplot(vlbob, ts, docur, opts, umax, vmax, &modified);
/*
 * If the data were modified, mark the map and beam as out of date.
 */
  if(vlbmap && modified)
    vlbmap->domap = vlbmap->dobeam = MAP_IS_STALE;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Plot visibility time sampling of telescopes of one or more sub-arrays.
 *
 * Input:
 *  isub   int   A single sub-array to be plotted.
 *  cif    int   The index of the start IF, or -1 for the default.
 */
static Template(timpl_fn)
{
  Subspec *ss=NULL; /* The specification of the initial sub-array */
  int cif = -1;     /* The default start IF */
  int modified=0;   /* True if data are edited in tplot */
  int ierr=0;       /* Error status returned by tplot. */
  char *opts="";    /* Default plot options */
/*
 * Make sure that we have data to be plotted.
 */
  if(nodata("tplot", OB_SELECT))
    return -1;
/*
 * Ensure that PGPLOT is open and clear the current display.
 */
  if(make_open() == -1)
    return -1;
/*
 * Arguments.
 */
  switch(npar) {
  case 2:
    cif = *INTPTR(invals[1]) - 1;
  case 1:
    ss = read_Subspec(vlbob, *STRPTR(invals[0]), NULL, 0);
    if(!ss)
      return -1;
  };
/*
 * Does the tflags options string contain any entries?
 * If so override the default options string.
 */
  if(**STRPTR(&tflags)=='\0')
    lprintf(stdout, "Using default options string \"%s\"\n", opts);
  else {
    opts = *STRPTR(&tflags);
    lprintf(stdout,
       "Overriding default options with user defined tflags=\"%s\"\n", opts);
  };
/*
 * Plot the data.
 */
  ierr = timplt(vlbob, ss, cif, 1, opts, &modified);
/*
 * Newly flagged or unflagged visibilities affect the map and beam.
 */
  if(vlbmap && modified)
    vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
  return ierr ? -1 : no_error;
}

/*.......................................................................
 * Determine and apply residual baseline corrections, one amp+phase
 * correction per baseline over the observation.
 * 
 * Input:
 *  basename   char *  Baseline specification string.
 */
static Template(resof_fn)
{
  Basespec *bs=NULL; /* Baseline specification descriptor */
  int cif;           /* The IF being corrected */
/*
 * Make sure that we have data to be corrected.
 */
  if(nodata("resoff", OB_SELECT))
    return -1;
/*
 * Check the baseline name given.
 */
  bs = read_Basespec(vlbob, npar>0 ? *STRPTR(invals[0]) : "", NULL, 0);
  if(!bs || next_base(vlbob, FIND_FIRST, 1, bs->nfix, 1, 0, 1, bs))
    return -1;
/*
 * Do we have a model to use?
 */
  if(vlbob->model==NULL) {
    lprintf(stderr, "resoff: No model to use in residual determination\n");
    return -1;
  };
/*
 * The change to the data calibration invalidates the current residual map.
 */
  if(vlbmap != NULL)
    vlbmap->domap = MAP_IS_STALE;
/*
 * Process all sampled IFs one at a time.
 */
  for(cif=0; (cif=nextIF(vlbob, cif, 1, 1)) >= 0; cif++) {
    Basespec bstmp = *bs;
    if(getIF(vlbob, cif))
      return -1;
    do {
      if(resoff(vlbob, 0, bstmp.base, bstmp.isub))
	return -1;
    } while(next_base(vlbob, FIND_NEXT, 1, bstmp.nfix, 1, 0, 0, &bstmp)==0);
  };
  return no_error;
}

/*.......................................................................
 * Remove all baseline based calibration offsets.
 */
static Template(unoff_fn)
{
  int dophs=0;  /* Default to no undoing of phase corrections */
  int doamp=0;  /* Default to no undoing of amplitude corrections */
/*
 * Make sure that there are data to be un-corrected.
 */
  if(nodata("clroff", OB_INDEX))
    return -1;
/*
 * Override defaults with user arguments (Note fallthrough).
 */
  switch(npar) {
  case 2:
    doamp = *LOGPTR(invals[1]);
  case 1:
    dophs = *LOGPTR(invals[0]);
  };
/*
 * NULL operation?
 */
  if(!doamp && !dophs)
    lprintf(stderr,"clroff: Neither amplitude nor phase specified by user\n");
  else {
    if(vlbmap)
      vlbmap->domap = MAP_IS_STALE;
/*
 * Clear the offsets.
 */
    if(clroff(vlbob, 1, doamp, dophs))
      return -1;
/*
 * Keep user informed about what has been done.
 */
    if(dophs)
      lprintf(stdout,
	     "clroff: All baseline phase corrections have been un-done.\n");
    if(doamp)
      lprintf(stdout,
	     "clroff: All baseline amplitude corrections have been un-done.\n");
  };
  return no_error;
}

/*.......................................................................
 * Save data, windows, model and restored map using standard extension
 * names.
 *
 * Input:
 *  prefix  char *  The name of the file to append extensions to.
 */
static Template(save_fn)
{
  static Descriptor filearg={'c',0,R_ONLY,1,{1,1,1},NULL};
  static Descriptor boolarg={'c',0,R_ONLY,1,{1,1,1},NULL};
  Descriptor *dsc[] = {&filearg, &boolarg};
  char *fname; /* Pointer to dynamic array to compose file names in */
  char *bname; /* Pointer to input base-name string */
  int slen;    /* The length of the input string */
  int hasmod;  /* True if at least one of the models contains compoenents */
  int ierr=0;  /* Error status flag */
/*
 * Check that there is something to be saved.
 */
  if(nodata("save", OB_INDEX))
    return -1;
/*
 * Get the pointer of the input base-name string and find its length.
 */
  bname = *STRPTR(invals[0]);
  slen = strlen(bname);
/*
 * Allocate a string to hold MAXSUF characters more than the length of
 * the base name.
 */
  fname = (char *) malloc((size_t) slen+MAXSUF);
  if(fname==NULL) {
    lprintf(stderr, "save: Insufficient memory to compose file names in\n");
    return -1;
  };
/*
 * Make the file name descriptor point to this string.
 */
  VOIDPTR(&filearg) = &fname;
/*
 * Copy the base name into the first slen elements.
 */
  strncpy(fname, bname, slen);
/*
 * Write the UV FITS file.
 */
  strcpy(&fname[slen], uvf_nam);
  ierr = wobs_fn(dsc, 1, outvals);
/*
 * Write the model file.
 */
  hasmod = vlbob->model->ncmp + vlbob->newmod->ncmp > 0;
  if(!ierr && hasmod) {
    strcpy(&fname[slen], mod_nam);
    ierr = wmodel_fn(dsc, 1, outvals);
  };
/*
 * Write the continuum-model file.
 */
  if(!ierr && vlbob->cmodel->ncmp + vlbob->cnewmod->ncmp > 0) {
    char docont=1;
    strcpy(&fname[slen], cmod_nam);
    VOIDPTR(&boolarg) = &docont;
    ierr = wmodel_fn(dsc, 2, outvals);
  };
/*
 * Write the windows file.
 */
  if(!ierr && vlbwins != NULL && vlbwins->nwin!=0) {
    strcpy(&fname[slen], win_nam);
    ierr = wwins_fn(dsc, 1, outvals);
  };
/*
 * Write a fits file of the restored map.
 */
  if(!ierr && vlbmap != NULL && hasmod) {
    strcpy(&fname[slen], fits_nam);
    ierr = wmap_fn(dsc, 1, outvals);
  };
/*
 * Write the table of selection-specific models.
 */
  if(!ierr && num_ModelTable_entries(vlbob->mtab)) {
    strcpy(&fname[slen], mtab_nam);
    ierr = write_models_fn(dsc, 1, outvals);
  };
/*
 * Write a command file that may subsequently be invoked to restore all
 * the files and running parameters to their current state.
 */
  if(!ierr) {
    strcpy(&fname[slen], par_nam);
    ierr = wrtpars(fname, bname);
  };
/*
 * Clean-up and return.
 */
  free(fname);
  if(ierr)
    return ierr;
  return no_error;
}

/*.......................................................................
 * Restore data, windows, and model using standard extension names.
 *
 * Input:
 *  prefix  char *  The name of the file to append extensions to.
 */
static Template(get_fn)
{
  static Descriptor filearg={'c',0,R_ONLY,1,{1,1,1},NULL};
  static Descriptor boolarg={'c',0,R_ONLY,1,{1,1,1},NULL};
  Descriptor *dsc[] = {&filearg, &boolarg};
  char *fname; /* Pointer to dynamic array to compose file names in */
  char *bname; /* Pointer to input base-name string */
  int slen;    /* The length of the input string */
  int ierr=0;  /* Error status flag */
/*
 * Get the pointer of the input base-name string and find its length.
 */
  bname = *STRPTR(invals[0]);
  slen = strlen(bname);
/*
 * Allocate a string to hold MAXSUF characters more than the length of
 * the base name.
 */
  fname = (char *) malloc((size_t) slen+MAXSUF);
  if(fname==NULL) {
    lprintf(stderr, "get: Insufficient memory to compose file names in\n");
    return -1;
  };
/*
 * Make the file name descriptor point to this string.
 */
  VOIDPTR(&filearg) = &fname;
/*
 * Copy the base name into the first slen elements.
 */
  strncpy(fname, bname, slen);
/*
 * Read the merge file.
 */
  strcpy(&fname[slen], uvf_nam);
  ierr = newob_fn(dsc, 1, outvals);
/*
 * Read the model file if it exist.
 */
  if(!ierr) {
    strcpy(&fname[slen], mod_nam);
    if(file_exists(fname)) {          /* Does the model file exist */
      ierr = rmodel_fn(dsc, 1, outvals);
    } else {
      lprintf(stdout, "Model file \"%s\" not available\n", fname);
    };
  };
/*
 * Read the continuum-model file if it exists.
 */
  if(!ierr) {
    strcpy(&fname[slen], cmod_nam);
    if(file_exists(fname)) {          /* Does the model file exist */
      char docont=1;
      VOIDPTR(&boolarg) = &docont;
      ierr = rmodel_fn(dsc, 2, outvals);
    };
  };
/*
 * Read the windows file.
 */
  if(!ierr) {
    strcpy(&fname[slen], win_nam);
    if(file_exists(fname)) {          /* Does the window file exist */
      ierr = rwins_fn(dsc, 1, outvals);
    } else {
      lprintf(stdout, "Window file \"%s\" not available\n", fname);
    };
  };
/*
 * Read the table of selection-specific models.
 */
  if(!ierr) {
    strcpy(&fname[slen], mtab_nam);
    if(file_exists(fname)) {          /* Does the multi-model file exist? */
      ierr = read_models_fn(dsc, 1, outvals);
    } else {
      lprintf(stdout, "Multi-model file \"%s\" not available\n", fname);
    };
  };
/*
 * Clean-up and return.
 */
  free(fname);
  if(ierr)
    return ierr;
  return no_error;
}

/*.......................................................................
 * Return an array of logarithmic contour levels.
 *
 * Input:
 *  absmin   float    The minimum +ve contour level.
 *  absmax   float    The maximum +ve contour level - default=100.0.
 *  factor   float    The exponential factor to use - default=2.
 * Output:
 *  mb_levs  float *  File-scope levels array.
 */
static Template(loglev_fn)
{
  const double tiny=1.0e-5;  /* Smallest allowable absmin,absmax or factor */
  const int maxlev=MAXARG*10;/* Largest allowable number of levels */
  int nlev;             /* The number of levels in the output array */
  float *levs;          /* Pointer into return array */
  double absmin;        /* Minimum contour level */
  double absmax=100.0;  /* Maximum contour level */
  double factor=2.0;    /* Exponential factor */
  double dtmp;
  int i;
/*
 * Get the arguments.
 */
  switch(npar) {
  case 3:
    factor = *FLTPTR(invals[2]);
  case 2:
    absmax = *FLTPTR(invals[1]);
  case 1:
    absmin = *FLTPTR(invals[0]);
    break;
  default:
    lprintf(stderr, "loglevs: Syserr - npar=%d\n", npar);
    return -1;
  };
/*
 * Ensure positivity.
 */
  if(absmin < 0.0f)
    absmin = -absmin;
  if(absmax < 0.0f)
    absmax = -absmax;
  if(factor < 0.0f)
    factor = -factor;
/*
 * Ensure that absmin < absmax.
 */
  if(absmin>absmax) {dtmp=absmin; absmin=absmax; absmax=dtmp;};
/*
 * Check arguments.
 */
  if(absmin<tiny || absmax<tiny) {
    lprintf(stderr, "loglevs: Bad limits min=%g max=%g\n", absmin, absmax);
    return -1;
  };
  if(factor<tiny || factor-1.0f < tiny) {
    lprintf(stderr, "loglevs: Illegal factor=%g\n", factor);
    return -1;
  };
/*
 * Calculate the number of levels in the output array.
 */
  nlev = 2 + log(absmax/absmin)/log(factor);
  if(nlev>maxlev) {
    lprintf(stderr, "loglevs: Too many levels (%d) required - maximum is %d\n",
	    nlev, maxlev);
    return -1;
  };
/*
 * Realloc the levels array for the new number of elements.
 */
  if(mb_levs.num_el < nlev) {
    valof_free(&mb_levs);
    mb_levs.num_el = mb_levs.adim[0] = 0;
    if( (VOIDPTR(&mb_levs) = valof_alloc(nlev, 'f')) == NULL)
      return -1;
    mb_levs.num_el = nlev;
  };
  mb_levs.adim[0] = nlev;
/*
 * Get a pointer to the return array.
 */
  levs = FLTPTR(&mb_levs);
/*
 * Make the first level a -ve version of the lower level limit.
 */
  levs[0] = -absmin;
/*
 * Build the exponential ramp of levels.
 */
  for(i=1; i<nlev; i++)
    levs[i] = absmin * pow(factor, (double) (i-1));
/*
 * Inform the user.
 */
  lprintf(stdout, "The new contour levels are:\n");
  for(i=0; i<nlev; i++)
    lprintf(stdout, " %g", levs[i]);
  lprintf(stdout, "\n");
  return no_error;
}

/*.......................................................................
 * Plot a the residual/clean map or beam and allow manipulation of
 * windows.
 */
static Template(maplot_fn)
{
  float *levs;      /* Pointer to array of contour levels */
  int nlevs;        /* The number of contours in the 'levs' array */
  int docont;       /* True if contours are required */
  int domap=1;      /* Default to plotting the map */
  int domod=0;      /* Default to not plotting the model components */
  int docln=0;      /* If 1 plot the restored map */
  int dovect=0;     /* True to plot a polarization map */

/* List valid argument values for the map-type argument */

  enum {PL_BEAM, PL_MAP, PL_CLEAN, PL_PMAP, PL_PCLN};
  static Enumpar images[] = {
    {"map",  PL_MAP},                       /* Plot the residual map */
    {"beam", PL_BEAM},                      /* Plot the beam */
    {"cln",  PL_CLEAN}, {"clean", PL_CLEAN},/* Plot the clean map */
    {"pmap", PL_PMAP},                      /* Plot residual polarization map */
    {"pcln", PL_PCLN},                      /* Plot clean polarization map */
  };
  static Enumtab *imtab=NULL; /* Symbol table of member enumerators */
/*
 * Check that a map and beam have been allocated and that an observation
 * has been read.
 */
  if(nomap("mapplot") || nodata("mapplot", OB_SELECT))
    return -1;
/*
 * Construct the enumerator symbol table if not already done.
 */
  if(!imtab && !(imtab=new_Enumtab(images, sizeof(images) / sizeof(Enumpar),
				   "mapplot image")))
    return -1;
/*
 * See if the map or beam descriptor was sent as the 1 optional argument.
 */
  if(npar > 0) {
/*
 * Lookup the image to be displayed.
 */
    Enumpar *image = find_enum(imtab, *STRPTR(invals[0]));
    if(!image)
      return -1;
/*
 * Arrange for the requested image to be displayed.
 */
    switch(image->id) {
    case PL_MAP:
      domap = 1;
      break;
    case PL_BEAM:
      domap = 0;
      break;
    case PL_CLEAN:
      domap = 1;
      docln = 1;
      break;
    case PL_PMAP:
      domap = 1;
      dovect = 1;
      break;
    case PL_PCLN:
      domap = 1;
      docln = 1;
      dovect = 1;
      break;
    };
  };
/*
 * See if the model components are to be plotted symbolically on the plot.
 */
  if(npar>1 && *LOGPTR(invals[1]))
    domod=1;
/*
 * Are polarization vectors are to be plotted?
 */
  if(dovect) {
/*
 * Check that polarization vectors have been enabled.
 */
    if(mappar.vect.scale == 0.0) {
      lprintf(stderr,
	      "Please use the 'polvec' command to configure the vectors.\n");
      return -1;
    };
/*
 * Make the polarization angle and intensity maps, placing them before
 * and after the main map in the vlbmap->map array.
 */
    if(make_polmap(docln) == -1)
      return -1;
/*
 * If a restored map is required, restore it unless it has already been
 * restored and not invalidated by other commands since restoring.
 * restore_fn uniquely marks the map as invalidating the
 * residual map by setting vlbmap->domap=RESTORED while all other
 * functions mark vlbmap->domap=1.
 */
  } else if(docln) {
    if(!vlbmap->ncmp || vlbmap->domap!=MAP_IS_CLEAN)
      if(restore_fn(invals,0,outvals) == -1)
	return -1;
  }
/*
 * If plotting residual map or dirty beam, re-invert if necessary.
 */
  else if( ((domap && vlbmap->domap) || (vlbmap->dobeam && !domap)) &&
	  invert_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Get the contour levels array.
 */
  if(mb_levs.adim[0] > 1) {
    levs = FLTPTR(&mb_levs);
    nlevs = mb_levs.adim[0];
  } else {
    levs = NULL;
    nlevs = 0;
  };
/*
 * If no windows have yet been set, create a new Mapwin container.
 */
  if(vlbwins == NULL)
    vlbwins = new_Mapwin();
/*
 * Ensure that PGPLOT is open.
 */
  if(make_open() == -1)
    return -1;
/*
 * Only allow contours on restored maps.
 */
  docont = mappar.docont && ((vlbmap->ncmp && domap) ||
			     mappar.ctab->cmap->class==CM_NONE);
/*
 * Plot the map or beam.
 */
  if(maplot(vlbob, vlbmap, vlbwins, &mappar.mpb, &mappar.vect,
	    domap, mappar.ctab, docont, dovect, domod, levs, nlevs,
	    mappar.cmul, mappar.box, mapmarkers))
    return -1;
  return no_error;
}

/*.......................................................................
 * Set/display limits on the displayed area in maplot.
 *
 * Input:
 *  xa   float  An upper or lower limit to the displayed X-extent.
 *  xb   float  The opposite limit to xa.
 *  ya   float  An upper or lower limit to the displayed Y-extent.
 *  yb   float  The opposite limit to ya.
 */
static Template(xyrange_fn)
{
/*
 * Install up to four arguments.
 */
  switch(npar) {
  case 4:
    mappar.box[3] = xytorad(*FLTPTR(invals[3]));
  case 3:
    mappar.box[2] = xytorad(*FLTPTR(invals[2]));
  case 2:
    mappar.box[1] = xytorad(*FLTPTR(invals[1]));
  case 1:
    mappar.box[0] = xytorad(*FLTPTR(invals[0]));
  };
/*
 * Report on the new limits.
 */
  lprintf(stdout, "The mapplot display area will be limited to:\n");
/*
 * The Right ascension plot limits.
 */
  lprintf(stdout, " Relative  RA: ");
  if(mappar.box[0]==mappar.box[1])
    lprintf(stdout, "(The whole available range)\n");
  else
    lprintf(stdout, "%.3g -> %.3g (%s)\n",
	    radtoxy(mappar.box[0]), radtoxy(mappar.box[1]), mapunits(U_TLAB));
/*
 * The declination plot limits.
 */
  lprintf(stdout, " Relative DEC: ");
  if(mappar.box[2]==mappar.box[3])
    lprintf(stdout, "(The whole available range)\n");
  else
    lprintf(stdout, "%.3g -> %.3g (%s)\n",
	    radtoxy(mappar.box[2]), radtoxy(mappar.box[3]), mapunits(U_TLAB));
  return no_error;
}

/*.......................................................................
 * Create a command file that contains the calls required to re-instate
 * the current operating parameters and files at a later date.
 *
 * Input:
 *  parname   char *   The name for the new command file.
 *  basename  char *   The base part of the parname.
 */
static int wrtpars(char *parname, char *basename)
{
  FILE *fp;      /* File descriptor of the command file */
  int waserr=0;  /* True after an error occurs */
/*
 * Parameters are meaningless without an observation.
 */
  if(nodata("wrtpars", OB_INDEX))
    return -1;
/*
 * Attempt to open the new command file.
 */
  if(parname) {
    fp = fopen(parname, "w");
    if(fp==NULL) {
      lprintf(stderr, "wrtpars: Error opening command file: %s\n", parname);
      return -1;
    };
    lprintf(stdout, "Writing difmap environment to: %s\n", parname);
/*
 * Write a time stamp for future reference.
 */
    waserr = waserr || lprintf(fp,
	     "! Command file created by the difmap on %s\n", date_str()) < 0;
  } else {
    fp = stdout;  /* Write to stdout if no file name was given */
  };
/*
 * Before doing anything else, set up the units with which map XY and UV
 * coordinates are to be refered by.
 */
  waserr = waserr || lprintf(fp, "mapunits %s\n", mapunits(U_NAME)) < 0;
/*
 * Arrange to enable or disable multi-model mode before any models
 * are read in.
 */
  waserr = waserr || lprintf(fp, "multi_model %s\n", multi_model_mode ?
			     "true" : "false") < 0;
/*
 * Store the "get" command to handle reading whatever data files are
 * available.
 */
  if(basename)
    waserr = waserr || lprintf(fp, "get %s\n", basename) < 0;
/*
 * Write the commands.
 */
  if(vlbmap) {
    waserr = waserr || lprintf(fp, "mapsize %d,%g, %d,%g\n",
		       vlbmap->nx, radtoxy(vlbmap->xinc),
		       vlbmap->ny, radtoxy(vlbmap->yinc)) < 0;
  };
/*
 * Write the select command, including all channel ranges currently selected.
 */
  if(ob_ready(vlbob, OB_SELECT, NULL)) {
    Chlist *cl = vlbob->stream.cl;
    int ir;
    waserr = waserr || lprintf(fp, "select %s",
			       Stokes_name(vlbob->stream.pol.type)) < 0;
    for(ir=0; ir<cl->nrange; ir++) {
      waserr = waserr || lprintf(fp, ", %d, %d",
				 cl->range[ir].ca+1, cl->range[ir].cb+1) < 0;
    };
    waserr = waserr || lprintf(fp, "\n") < 0;
  };
  waserr = waserr || lprintf(fp, "uvtaper %g, %g\n",
			     invpar.gauval, wavtouv(invpar.gaurad)) < 0;
  waserr = waserr || lprintf(fp, "uvrange %g, %g\n",
			     wavtouv(invpar.uvmin), wavtouv(invpar.uvmax)) < 0;
  waserr = waserr || lprintf(fp, "uvweight %g, %g, %s\n",
			     invpar.uvbin, invpar.errpow,
			     invpar.dorad?"true":"false") < 0;
  if(vlbob->uvzero.wt > 0.0f) {
    waserr = waserr || lprintf(fp, "uvzero %g, %g\n", vlbob->uvzero.amp,
			       vlbob->uvzero.wt) < 0;
  };
  waserr = waserr || lprintf(fp, "selftaper %g, %g\n",
			     slfpar.gauval, wavtouv(slfpar.gaurad)) < 0;
  waserr = waserr || lprintf(fp, "selflims %g, %g\n",
			     slfpar.maxamp, slfpar.maxphs * rtod) < 0;
  waserr = waserr || lprintf(fp, "xyrange %g, %g, %g, %g\n",
			     radtoxy(mappar.box[0]), radtoxy(mappar.box[1]),
			     radtoxy(mappar.box[2]), radtoxy(mappar.box[3]))<0;
  waserr = waserr || lprintf(fp, "beamloc %g, %g, %g, %g\n",
			     mappar.mpb.xc, mappar.mpb.yc,
			     mappar.mpb.minsize, mappar.mpb.maxsize)<0;
  waserr = waserr || lprintf(fp, "polvec %g, %g, %g, %d, %d\n",
			     radtoxy(mappar.vect.scale), mappar.vect.icut,
			     mappar.vect.pcut, mappar.vect.dx,
			     mappar.vect.dy) < 0;
  waserr = waserr || lprintf(fp, "integer niter; niter=%d\n", clnpar.niter)<0;
  waserr = waserr || lprintf(fp, "float gain; gain=%g\n", clnpar.gain)<0;
  waserr = waserr || lprintf(fp, "float cutoff; cutoff=%g\n", clnpar.cutoff)<0;
  waserr = waserr || lprintf(fp, "float cmul; cmul=%g\n", mappar.cmul)<0;
  waserr = waserr || lprintf(fp, "logical docont; docont=%s\n",
			     mappar.docont ? "true" : "false") < 0;
  {
    Ctable *ctab = mappar.ctab;
    waserr = waserr || lprintf(fp, "mapcolor %s, %f, %f\n",
			       ctab->cmap->class==CM_GREY ? "grey":"color",
			       ctab->contra, ctab->bright) < 0;
    waserr = waserr || lprintf(fp, "mapfunc %s, %f, %f\n",
		       name_Cmtran(ctab->tran), ctab->vmin, ctab->vmax) < 0;
  };
  waserr = waserr || lprintf(fp, "string vflags; vflags=\"%s\"\n", *STRPTR(&vflags)) < 0;
  waserr = waserr || lprintf(fp, "string rflags; rflags=\"%s\"\n", *STRPTR(&rflags)) < 0;
  waserr = waserr || lprintf(fp, "string pflags; pflags=\"%s\"\n", *STRPTR(&pflags)) < 0;
  waserr = waserr || lprintf(fp, "string tflags; tflags=\"%s\"\n", *STRPTR(&tflags)) < 0;
  waserr = waserr || lprintf(fp, "selfflag %s, %d, %d\n",
			     slfpar.doflag ? "true":"false",
			     slfpar.p_mintel, slfpar.a_mintel) < 0;
  if(vlbob->geom.east != 0.0 || vlbob->geom.north != 0.0) {
    waserr = waserr || lprintf(fp, "shift %g, %g\n",
                    radtoxy(vlbob->geom.east), radtoxy(vlbob->geom.north)) < 0;
  };
/*
 * Write commands to record non-unity selfcal antenna weights and
 * antenna constraint flags.
 */
  if(!waserr) {
    Subarray *sub = vlbob->sub;
    int isub;
    for(isub=0; !waserr && isub<vlbob->nsub; isub++,sub++) {
      Station *tel = sub->tel;
      int itel;
      for(itel=0; !waserr && itel<sub->nstat; itel++,tel++) {
	if(fabs(tel->antwt-1.0f) > 0.01  ||  tel->antfix) {
	  waserr = lprintf(fp, "selfant \"%d:%s\", %s, %g\n", isub+1,
		       tel->name, tel->antfix?"true":"false", tel->antwt) < 0;
	};
      };
    };
  };
/*
 * Write commands to record the interscan gap of each sub-array.
 */
  if(!waserr) {
    Subarray *sub;
    int isub;
    int same=1;  /* True if all interscan gaps are identical */
/*
 * Determine whether all the interscan gaps are identical.
 */
    sub = vlbob->sub;
    for(isub=0; !waserr && same && isub<vlbob->nsub; isub++)
      same = sub[isub].scangap == sub[0].scangap;
/*
 * If all gaps are identical then only one statement is required.
 */
    if(same) {
      waserr = lprintf(fp, "scangap %g\n", sub[0].scangap) < 0;
    } else {
      for(isub=0; !waserr && isub<vlbob->nsub; isub++)
	waserr = lprintf(fp, "scangap %g, %d\n", sub[isub].scangap, isub+1) < 0;
    };
  };
/*
 * Write out float array assignments.
 */
  waserr = waserr || w_flt_array(fp, "levs", &mb_levs);
/*
 * Write the commands that are needed to restore the list of map markers.
 */
  waserr = waserr || write_marker_commands(fp);
/*
 * Report errors.
 */
  if(parname && waserr)
    lprintf(stderr, "wrtpars: Error writing parameters to: %s\n", parname);
/*
 * Close the command file.
 */
  if(fp!=stdout)
    fclose(fp);
  return waserr ? -1 : no_error;
}

/*.......................................................................
 * Private function used by wrtpars() to write a float array assignment.
 *
 * Input:
 *  fp                  The stream to write to.
 *  name       char *   The name of the variable.
 *  dsc  Descriptor *   Variable descriptor.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int w_flt_array(FILE *fp, char *name, Descriptor *dsc)
{
  int waserr=0;               /* True after I/O error */
  int ntotal = dsc->adim[0];  /* Number of elements to write */
  float *array = FLTPTR(dsc); /* Pointer into float array */
  int nnew;                   /* The number to be written in one assignement */
  int start;                  /* The First index of an assignment */
  int i;
/*
 * Ensure that the array has been declared to take the number of elements
 * that are to be assigned.
 */
  waserr = lprintf(fp, "float %s(%d)\n", name, ntotal) < 0;
/*
 * Assign up to MAXARG elements at a time.
 */
  for(start=0; !waserr && start<ntotal; start += MAXARG) {
/*
 * How many elements to be assigned in this statement?
 */
    nnew = ntotal-start;
    if(nnew >= MAXARG)
      nnew = MAXARG;
/*
 * Write the assignment statement.
 */
    waserr = lprintf(fp, "%s(%d:%d) =", name, start+1, start+nnew) < 0;
/*
 * Write the values to be assigned - up to MAXARG elements at a time.
 */
    for(i=0; !waserr && i<nnew; i++)
      waserr = lprintf(fp, " %g%c", *array++, (i<nnew-1 ? ',':'\n')) < 0;
  };
  return waserr;
}

/*.......................................................................
 * Control the identification and treatment of un-correctable visibilities
 * in sub-sequent self-calibration commands.
 *
 * Input:
 *  doflag    char  If true flag un-correctable visibilities.
 *  p_mintel   int  The minimum number of telescopes required for phase
 *                  solutions.
 *  a_mintel   int  The minimum number of telescopes required for amplitude
 *                  solutions.
 */
static Template(sflag_fn)
{
/*
 * The flagging is undefined until an observation has been read.
 */
  if(nodata("selfflag", OB_INDEX))
    return -1;
/*
 * Override existing values where given.
 */
  if(npar > 0) {
    switch (npar) {
    case 3:
      slfpar.a_mintel = *INTPTR(invals[2]);
    case 2:
      slfpar.p_mintel = *INTPTR(invals[1]);
    case 1:
      slfpar.doflag = *LOGPTR(invals[0]);
    };
/*
 * Enforce limits - it is impossible to have closed arrays containing
 * less than 3 telescopes. In such cases substitute 0 to turn off checks.
 */
    if(slfpar.a_mintel < 3)
      slfpar.a_mintel = 0;
    if(slfpar.p_mintel < 3)
      slfpar.p_mintel = 0;
  };
/*
 * Report disposition of effected visibilities in phase-only self-cal.
 */
  if(slfpar.p_mintel > 0) {
    lprintf(stdout, "- In phase-only self-cal, good data on baselines that are not in closed\n");
    lprintf(stdout, "  arrays of at least %d telescopes will %s.\n",
	    slfpar.p_mintel, (slfpar.doflag ? "be flagged" : "not be used"));
  } else {
    lprintf(stdout, "- In phase-only self-cal, all un-flagged data will be used\n");
  };
/*
 * Report disposition of effected visibilities in amplitude self-cal.
 */
  if(slfpar.a_mintel > 0) {
    lprintf(stdout, "- In amplitude self-cal, good data on baselines that are not in closed\n");
    lprintf(stdout, "  arrays of at least %d telescopes will %s.\n",
	    slfpar.a_mintel, (slfpar.doflag ? "be flagged" : "not be used"));
  } else {
    lprintf(stdout, "- In amplitude self-cal, all un-flagged data will be used\n");
  };
  return no_error;
}

/*.......................................................................
 * Set and/or display extra telescope based weights and flags, to be
 * used subsequently in slfcal().
 *
 * Input:
 *  tname     char *  The name of the station to be addressed. If no name
 *                    is provided, all telescopes are addressed.
 *  fix       char    Logical argument. If true the gain of the addressed
 *                    station(s) will not be allowed to vary.
 *  weight   float    The extra weight multiplier to apply to the
 *                    addressed telescope(s).
 */
static Template(selfant_fn)
{
  Telspec *ts=NULL;  /* Telescope specification descriptor */
  float weight=0.0f; /* The telescope weight to be applied */
  int dofix=0;       /* If true fix the gain of the addressed telescope */
  char *tname="";    /* The telescope name */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("selfant", OB_INDEX))
    return -1;
/*
 * Parse arguments - fall-through intentional.
 */
  switch(npar) {
  case 3:
    weight = *FLTPTR(invals[2]);
    if(weight <= 0.0f) {
      lprintf(stderr, "selfant: Weights must be greater than zero.\n");
      return -1;
    };
  case 2:
    dofix = *LOGPTR(invals[1]);
  case 1:
    tname = *STRPTR(invals[0]);
  };
/*
 * Decode the telescope specification.
 */
  ts = read_Telspec(vlbob, tname, NULL, 0);
  if(!ts || next_tel(vlbob, FIND_FIRST, 1, ts->nfix, 0, 1, ts))
    return -1;
/*
 * Iterate through the specified telescope(s).
 */
  do {
    Station *tel = vlbob->sub[ts->isub].tel + ts->ta;
/*
 * Apply any selections that were provided.
 */
    if(npar > 1) {
      tel->antfix = dofix;
      if(weight>0.0f)
	tel->antwt = weight;
    };
/*
 * Report the selections.
 */
    lprintf(stdout, "%d:%-8s self-cal  status=%s  weight=%g\n",
	    ts->isub+1, tel->name, tel->antfix ? "fixed" : "correctable",
	    tel->antwt);
  } while(next_tel(vlbob, FIND_NEXT, 1, ts->nfix, 0, 0, ts)==0);
  return no_error;
}

/*.......................................................................
 * Display the history lines of the current observation.
 */
static Template(hist_fn)
{
/*
 * Make sure that an observation has been read.
 */
  if(nodata("showhist", OB_INDEX))
    return -1;
/*
 * Invoke the history displayer.
 */
  return showhist(vlbob, 1) ? -1 : no_error;
}

/*.......................................................................
 * Allow selection of a new UV stream.
 *
 * Input:
 *  stokes  char *  The name of a known stokes parameter.
 *  bchan    int    Index of the start channel of the first channel range.
 *  echan    int    Index of the last channel of the first channel wanted.
 *  bchan,echan...  Further channel ranges.
 */
static Template(uvsel_fn)
{
  Stokes stokes=NO_POL; /* Enumerated stokes type */
  Chlist *cl = NULL;    /* NULL channel list re-selects original list */
  int bchan = 0;        /* 0-relative index of first channel */
  int echan = 0;        /* 0-relative index of last channel */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("select", OB_INDEX))
    return -1;
/*
 * Set defaults equal to the currently selected stream parameters.
 */
  stokes = vlbob->stream.pol.type;
/*
 * Get the stokes argument.
 */
  if(npar > 0) {
    stokes = Stokes_id(*STRPTR(invals[0]));
    if(stokes==NO_POL)
      return -1;
  };
/*
 * Get channel ranges if given.
 */
  if(npar > 1) {
    int par;
/*
 * Allocate a channel-range container in which to record each channel
 * range.
 */
    if((cl = new_Chlist())==NULL)
      return -1;
/*
 * Interpret paired channel-range arguments.
 */
    for(par=1; par<npar; par += 2) {
      bchan = *INTPTR(invals[par]) - 1;
/*
 * If the last channel range argument is unpaired, assume that it describes
 * a single channel and set echan = bchan.
 */
      echan = par+1 < npar ? (*INTPTR(invals[par+1]) - 1) : bchan;
      if(add_crange(cl, bchan, echan)) {
	del_Chlist(cl);
	return -1;
      };
    };
  };
/*
 * Mark the map as out of date.
 */
  if(vlbmap)
    vlbmap->domap = vlbmap->dobeam = MAP_IS_STALE;
/*
 * Instate the selected stream.
 */
  if(ob_select(vlbob, !multi_model_mode, cl, stokes))
    return -1;
  return no_error;
}

/*.......................................................................
 * This function should be called by all functions that use vlbob.
 * If vlbob is not ready for use it will emit an error message that
 * includes the name of the calling function, and returns 1. The
 * calling function should then perform an error return.
 *
 * For example:
 *  if(nodata("invert", OB_SELECT))
 *    return -1;
 *
 * Input:
 *  cname  const char *  The command name of the calling function.
 *  state     Obstate    The initialization state required of the
 *                       observation. See obs.h.
 * Output:
 *  return        int    0 - OK.
 *                       1 - vlbob is not ready for use.
 */
static int nodata(const char *cname, Obstate state)
{
  char *message = NULL;   /* An error message */
/*
 * Check whether an observation has been read.
 */
  if(vlbob == NULL) {
    message = "No UV data has been read in yet - use the 'observe' command";
  } else if(!ob_ready(vlbob, state, cname)) {
    switch(vlbob->state) {
    case OB_INDEX:
      message = "Use the 'select' command to select a data stream";
      break;
    default:
      message = "Corrupt observation discarded";
      vlbob = del_Observation(vlbob);
      break;
    };
  };
/*
 * Was an error detected?
 */
  if(message) {
    lprintf(stderr, "%s: %s.\n", cname, message);
    return 1;
  };
/*
 * The observation descriptor vlbob is apparently usable.
 */
  return 0;
}

/*.......................................................................
 * This function should be called by all functions that use vlbmap.
 * If vlbmap is not ready for use it will emit an error message that
 * includes the name of the calling function, and returns 1. The
 * calling function should then perform an error return.
 *
 * For example:
 *  if(nomap("mapplot"))
 *    return -1;
 *
 * Input:
 *  cname  const char *  The command name of the calling function.
 * Output:
 *  return        int    0 - OK.
 *                       1 - vlbob is not ready for use.
 */
static int nomap(const char *cname)
{
/*
 * Check whether a map/beam container has been allocated.
 */
  if(vlbmap == NULL) {
    lprintf(stderr,
        "%s: No map or beam yet created - use the 'mapsize' command.\n", cname);
    return 1;
  };
/*
 * The map/beam container is usable.
 */
  return 0;
}

/*.......................................................................
 * Optionally change the scale factor applied to the UV FITS visibility
 * weights.
 *
 * Input:
 *  scale   float   The new scale factor to apply.
 * Output:
 *  return  float   Optional return value gives the current scale factor.
 */
static Template(wtscal_fn)
{
/*
 * Do we have weights to scale?
 */
  if(nodata("wtscale", OB_INDEX))
    return -1;
/*
 * Get the arguments.
 */
  if(npar > 0) {
/*
 * Mark the map as out of date.
 */
    if(vlbmap)
      vlbmap->domap = vlbmap->dobeam = MAP_IS_STALE;
/*
 * Apply the given incremental weight factor.
 */
    if(wtscale(vlbob, *FLTPTR(invals[0])))
      return -1;
    lprintf(stdout, "Changed weight scale factor to: %g\n",vlbob->geom.wtscale);
  } else {
    lprintf(stdout, "Current weight scale factor is: %g\n",vlbob->geom.wtscale);
  };
/*
 * Return the current weight scale factor if a return value is required.
 */
  if(outvals)
    *FLTPTR(outvals) = vlbob->geom.wtscale;
  return no_error;
}

/*.......................................................................
 * Return the flux, and X and Y positions of the max absolute map
 * values.
 *
 * Input:
 *  member    char * The name of a pixel attribute (x,y or flux).
 *  type      char * The type of pixel wanted (min,max,abs).
 * Output:
 *  return   float   The X,Y or flux value requested.
 */
static Template(peak_fn)
{
/*
 * List valid argument values for the pixel type.
 */
  enum {MAXPIX, MINPIX, ABSPIX};
  static Enumpar modes[] = {{"max", MAXPIX}, {"min", MINPIX}, {"abs",ABSPIX}};
  static Enumtab *modetab=NULL; /* Symbol table of mode enumerators */
  char *modename = "abs";    /* Default mode argument */
  Enumpar *mode;             /* Pointer to element of modes[] */

/* List valid argument values for the pixel member */

  enum {XPIX, YPIX, FPIX, PIXRA, PIXDEC};
  static Enumpar members[] = {{"x", XPIX}, {"y", YPIX}, {"flux",FPIX},
			      {"ra", PIXRA}, {"dec", PIXDEC}};
  static Enumtab *memtab=NULL; /* Symbol table of member enumerators */
  char *membername = "flux";   /* Default member argument */
  Enumpar *member;             /* Pointer to element of members[] */
  Mappix *mpix;            /* Descriptor of absolute max valued pixel */
/*
 * Construct the enumerator symbol tables if not already done.
 */
  if(!modetab &&
    !(modetab=new_Enumtab(modes,sizeof(modes)/sizeof(Enumpar),
			  "peak: type")))
    return -1;
  if(!memtab &&
    !(memtab=new_Enumtab(members,sizeof(members)/sizeof(Enumpar),
			 "peak: attribute")))
    return -1;
/*
 * Do we have a map?
 */
  if(nomap("peak"))
    return -1;
/*
 * Invert the map if it is not ready for use?
 */
  if(vlbmap->domap && vlbmap->domap != MAP_IS_CLEAN &&
     invert_fn(invals,0,outvals) == -1)
    return -1;
/*
 * For the default to be optionally overriden below, get the descriptor
 * of the pixel with the max abolute value.
 */
  mpix = fabs(vlbmap->maxpix.value) > fabs(vlbmap->minpix.value) ?
    &vlbmap->maxpix : &vlbmap->minpix;
/*
 * Overide defaults with user arguments.
 */
  switch(npar) {
  case 2:
    modename = *STRPTR(invals[1]);
  case 1:
    membername = *STRPTR(invals[0]);
  };
/*
 * Determine which pixel-type is required.
 */
  mode = find_enum(modetab, modename);
  if(!mode)
    return -1;
/*
 * Get the pointer of the associated pixel descriptor.
 */
  switch(mode->id) {
  case MAXPIX:
    mpix = &vlbmap->maxpix;
    break;
  case MINPIX:
    mpix = &vlbmap->minpix;
    break;
  case ABSPIX: default:
    mpix = fabs(vlbmap->maxpix.value) > fabs(vlbmap->minpix.value) ?
      &vlbmap->maxpix : &vlbmap->minpix;
    break;
  };
/*
 * Determine the type of pixel member required.
 */
  member = find_enum(memtab, membername);
  if(!member)
    return -1;
/*
 * Assign the requested member as the function return value.
 */
  switch(member->id) {
  case XPIX:
    *FLTPTR(outvals) = radtoxy(mpix->xpos);
    break;
  case YPIX:
    *FLTPTR(outvals) = radtoxy(mpix->ypos);
    break;
  case FPIX: 
  default:
    *FLTPTR(outvals) = mpix->value;
    break;
  case PIXRA:
    *FLTPTR(outvals) = mpix->ra * rtod;
    break;
  case PIXDEC:
    *FLTPTR(outvals) = mpix->dec * rtod;
    break;
  };
/*
 * Done.
 */
  return no_error;
}

/*.......................................................................
 * Place a new CLEAN window around the max absolute flux in the current
 * map unless the position is already within a CLEAN window.
 *
 * Input:
 *  size  float  The size of the clean window wrt the fwhm aspect of
 *               the clean beam.
 *  doabs   int  If true search for the max absolute peak. Otherwise
 *               search out the max positive peak.
 */
static Template(pwin_fn)
{
  float size = 1.0f;  /* The relative size of the window */
  int doabs = 0;      /* By default search for the most positive peak */
/*
 * Has a map been alocated?
 */
  if(nomap("peakwin"))
    return -1;
/*
 * Get optional user arguments.
 */
  switch(npar) {
  case 2:
    doabs = *LOGPTR(invals[1]);
  case 1:
    size = *FLTPTR(invals[0]);
  };
/*
 * If the residual map or dirty beam are not ready for use, invert them.
 * The map is required to be up to date such that the peak pixel info
 * will be current. The beam is required to be up to date because the
 * estimated beam size from the last time that the beam was inverted
 * are used to set the window size.
 */
  if((vlbmap->domap || vlbmap->dobeam) && invert_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Create an empty Mapwin list if one hasn't been created yet.
 */
  if(vlbwins==NULL && (vlbwins=new_Mapwin()) == NULL)
    return -1;
/*
 * Ensure that the peak is windowed.
 */
  if(peakwin(vlbmap, vlbwins, size, doabs))
    return -1;
  return no_error;
}

/*.......................................................................
 * Fit components of the established and latest clean models that have
 * free parameters, to the visibility residuals.
 *
 * Input:
 *  niter   int   The max number of iterations.
 */
static Template(modfit_fn)
{
  int niter=0;   /* The max number of iterations */
/*
 * Anything to fit to?
 */
  if(nodata("modelfit", OB_SELECT))
    return -1;
/*
 * Get the number of iterations.
 */
  switch(npar) {
  case 1:
    niter = *INTPTR(invals[0]);
  };
/*
 * Check user arguments.
 */
  if(niter < 0) {
    lprintf(stderr, "modelfit: niter must be >= 0\n");
    return -1;
  };
/*
 * Mark the map as out of date.
 */
  if(vlbmap)
    vlbmap->domap = MAP_IS_STALE;
/*
 * Attempt to fit the variable part of the established and tentative models
 * to the residual visibilities.
 */
  if(fituvmodel(vlbob, niter, invpar.uvmin, invpar.uvmax))
    return -1;
  return no_error;
}

/*.......................................................................
 * Allow the user to edit selected parts of the tentative and established
 * models in an external editor.
 *
 * Input:
 *  dovar    int   If false just edit the variable part of the model.
 *                 If true, edit the full model.
 */
static Template(edmod_fn)
{
  int dovar = 0;  /* If false edit just the variable part of the models */
/*
 * No data to associate the model with?
 */
  if(nodata("edmod", OB_INDEX))
    return -1;
/*
 * Mark the map as out of date.
 */
  if(vlbmap)
    vlbmap->domap = MAP_IS_STALE;
/*
 * Edit the full model?
 */
  dovar = npar > 0 && *LOGPTR(invals[0]);
/*
 * Send the variable part of the model to an external editor.
 */
  if(obedmod(vlbob, dovar))
    return -1;
  return no_error;
}

/*.......................................................................
 * Display closure phases.
 *
 * Input:
 *  nrow     int   The number of plots per page.
 *  clsstr  char * The description of the triangle(s) to plot.
 *  cif      int   The start IF, or 0 for the default.
 *  npage    int   The max number of pages to plot when non-interactive.
 */
static Template(cpplt_fn)
{
  int nrow=0;        /* Default number of rows */
  int cif = -1;      /* Default start IF */
  int modified=0;    /* Returned true by clsplot() if the data were modified */
  int iret;          /* Return code from cpplot() */
  char *opts="efbm"; /* Default plot options */
  Trispec *ts=NULL;  /* Closure triangle designation */
  int npage = 0;     /* The max number of pages to plot when non-interactive. */
/*
 * Make sure that we have data to be plotted.
 */
  if(nodata("cpplot", OB_SELECT))
    return -1;
/*
 * Make sure that a PGPLOT device is open.
 */
  if(make_open() == -1)
    return -1;
/*
 * Check arguments.
 */
  switch (npar) {
  case 4:
    npage = *INTPTR(invals[3]);    
  case 3:
    cif = *INTPTR(invals[2]) - 1;
  case 2:
    ts = read_Trispec(vlbob, *STRPTR(invals[1]), NULL, 0);
    if(ts==NULL)
      return -1;
  case 1:
    nrow = *INTPTR(invals[0]);
  };
/*
 * Does the vflags options string contain any entries?
 * If so override the default options string.
 */
  if(**STRPTR(&vflags)=='\0')
    lprintf(stdout, "Using default options string \"%s\"\n", opts);
  else {
    opts = *STRPTR(&vflags);
    lprintf(stdout,
       "Overriding default options with user defined vflags=\"%s\"\n", opts);
  };
/*
 * Establish the latest clean model so that we have an up-to-date
 * representation of the UV model to display.
 */
  if(keep_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Plot the data.
 */
  iret = clsplot(vlbob, ts, cif, nrow, npage, 1, opts, 0, 0, 0, 0, &modified);
/*
 * Newly flagged or unflagged visibilities affect the map and beam.
 * If clsplot() reported the data as modified, mark both the map and beam
 * as out of date.
 */
  if(vlbmap && modified)
    vlbmap->dobeam = vlbmap->domap = MAP_IS_STALE;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Append a new line of history.
 *
 * Input:
 *  hline   char *  The history line to be written.
 */
static Template(addhis_fn)
{
  char *hline = *STRPTR(invals[0]);
/*
 * Do we have an observation to add history to?
 */
  if(nodata("addhist", OB_INDEX))
    return -1;
/*
 * Inform the user if the history line is over-length.
 */
  if(strlen(hline) > 80)
    lprintf(stdout, "History line truncated to 80 characters.\n");
/*
 * Add the new line of history.
 */
  if(add_hist(vlbob, hline))
    return -1;
  return no_error;
}

/*.......................................................................
 * Clear the currently recorded history.
 */
static Template(clrhis_fn)
{
  if(nodata("clrhist", OB_INDEX))
    return -1;
/*
 * Clear the currently recorded history.
 */
  if(clr_hist(vlbob))
    return -1;
  return no_error;
}

/*.......................................................................
 * Set the interscan separator gap for all sub-arrays, or a given sub-array.
 *
 * Input:
 *  gap   float  The gap to set (seconds). Zero sets the default.
 *  isub    int  The sub-array to assign the gap to. (0 selects all).
 */
static Template(scangap_fn)
{
  float gap = 0.0f;  /* The interscan gap */
  int isub = -1;     /* The sub-array to assign to (-1 -> all sub-arrays) */
/*
 * Check that we have data to assign the scan gaps to.
 */
  if(nodata("scangap", OB_INDEX))
    return -1;
/*
 * Are we changing anything?
 */
  if(npar > 0) {
/*
 * Override the argument defaults.
 */
    switch(npar) {
    case 2:
      isub = *INTPTR(invals[1]) - 1;
    case 1:
      gap = *FLTPTR(invals[0]);
    };
/*
 * Assign the new gap.
 */
    if(scangap(vlbob, gap, isub))
      return -1;
/*
 * If no arguments were specified, list the current values.
 */
  } else {
    Subarray *sub;
    int isub;
    int same=1;  /* True if all interscan gaps are identical */
/*
 * Determine whether all the interscan gaps are identical.
 */
    sub = vlbob->sub;
    for(isub=0; same && isub<vlbob->nsub; isub++)
      same = sub[isub].scangap == sub[0].scangap;
/*
 * If all gaps are identical then only one statement is required.
 */
    if(same) {
      lprintf(stdout,
         "The delimiting interscan gap is %g seconds in all sub-arrays.\n",
         sub[0].scangap);
    } else {
      for(isub=0; isub<vlbob->nsub; isub++)
	lprintf(stdout,
               "The delimiting interscan gap is %g seconds in sub-array %d.\n",
	       sub[isub].scangap, isub+1);
    };
  };
  return no_error;
}

/*.......................................................................
 * Set the units used to refer to map XY dimensions.
 *
 * Input:
 *  unit     char *   The official name of the unit.
 */
static Template(munit_fn)
{
/*
 * Set new map units?
 */
  if(npar > 0 && skyunits(*STRPTR(invals[0])))
     return -1;
/*
 * Report the currently selected map units.
 */
  lprintf(stdout,"Distances in the map plane now have units: %s.\n",
	  mapunits(U_TLAB));
  lprintf(stdout,"Distances in the UV plane now have units: %s.\n",
	  uvwunits(U_TLAB));
  return no_error;
}

/*.......................................................................
 * Add a new model component.
 *
 * Input:
 *  flux   float   Value of flux, plus whether variable.
 *  var     bool   True if the flux is a free parameter.
 *  x      float   Value of x, plus whether variable.
 *  y      float   Value of y, plus whether variable.
 *  var     bool   True if x and y are to be free parameters.
 *  major  float   Value of major, plus whether variable.
 *  var     bool   True if major is to be a free parameter.
 *  ratio  float   Value of ratio, plus whether variable.
 *  var     bool   True if ratio is to be a free parameter.
 *  phi    float   Value of phi, plus whether variable.
 *  var     bool   True if phi is to be a free parameter.
 *  type     int   The component type.
 *  freq0  float   The reference frequency (Hz).
 *  spcind float   The spectral index.
 *  var     bool   True if spcind is to be a free parameter.
 */
static Template(addmc_fn)
{
  int freepar=0;    /* Bitmap to specify which parameters are variable */
  float flux=0.0f;  /* The component flux */
  float x=0.0f;     /* The component Y-axis position (radians) */
  float y=0.0f;     /* The component X-axis position (radians) */
  float major=0.0f; /* The component major axis extent (radians) */
  float ratio=1.0f; /* The component axial ratio */
  float phi=0.0f;   /* The component major-axis position angle (radians) N->E */
  int type = 0;     /* Component type (default=delta). */
  float freq0=0.0f; /* The reference frequency */
  float spcind=0.0f;/* The spectral index */
/*
 * Check that we have an observation to add the model component to.
 */
  if(nodata("addcmp", OB_INDEX))
    return -1;
/*
 * Get user arguments.
 */
  switch(npar) {
  case 15:
    if(*LOGPTR(invals[14])) freepar |= M_SPCIND;
  case 14:
    spcind = *FLTPTR(invals[13]);
  case 13:
    freq0 = *FLTPTR(invals[12]);
  case 12:
    type = *INTPTR(invals[11]);
  case 11:
    if(*LOGPTR(invals[10])) freepar |= M_PHI;
  case 10:
    phi = *FLTPTR(invals[9]) * dtor;
  case 9:
    if(*LOGPTR(invals[8])) freepar |= M_RATIO;
  case 8:
    ratio = *FLTPTR(invals[7]);
  case 7:
    if(*LOGPTR(invals[6])) freepar |= M_MAJOR;
  case 6:
    major = xytorad(*FLTPTR(invals[5]));
  case 5:
    if(*LOGPTR(invals[4])) freepar |= M_CENT;
  case 4:
    y = xytorad(*FLTPTR(invals[3]));
  case 3:
    x = xytorad(*FLTPTR(invals[2]));
  case 2:
    if(*LOGPTR(invals[1])) freepar |= M_FLUX;
  case 1:
    flux = *FLTPTR(invals[0]);
    break;
  default:
    lprintf(stderr, "addcmp: Unexpected number of arguments: %d.\n", npar);
    return -1;
  };
/*
 * No type specified?
 * Substitute implicit type.
 */
  if(npar < 12)
    type = (npar < 6 || major == 0.0) ? M_DELT : M_GAUS;
/*
 * Check the component type.
 */
  if(type < 0 || type > 6) {
    lprintf(stderr, "addcmp: Unknown component type: %d.\n", type);
    return -1;
  };
/*
 * Add the component to the tentative model.
 */
  if(add_xycmp(vlbob->newmod, 1, freepar, flux, x, y, major, ratio, phi, type,
	       freq0, spcind) == NULL)
    return -1;
  return no_error;
}

/*.......................................................................
 * Return UV plane model and visibility statistics.
 *
 * Input:
 *  type     char *  The type of statistic to be returned.
 *                    rms    -  RMS deviation between model and data.
 *                    chisqr -  Chi squared.
 *                    nvis   -  Number of un-flagged visibilities.
 * Output:
 *  return  float    The requested statistic.
 */
static Template(uvstat_fn)
{
  enum {UVRMS, UVCHI, UVNVIS};
  static Enumpar types[] = {{"rms", UVRMS}, {"chisq", UVCHI}, {"nvis",UVNVIS}};
  static Enumtab *typtab = NULL;
  char *typename;        /* Statistic type name */
  Enumpar *type;         /* Pointer to element of types[] */
  Moddif md;
/*
 * Construct the enumerator symbol table if necessary.
 */
  if(!typtab &&
     !(typtab=new_Enumtab(types,sizeof(types)/sizeof(Enumpar),"UV statistic")))
    return -1;
/*
 * Do we have data?
 */
  if(nodata("uvstat", OB_SELECT))
    return -1;
/*
 * Get the user's choice of statistic.
 */
  typename = *STRPTR(invals[0]);
/*
 * Determine the type of statistic member required.
 */
  type = find_enum(typtab, typename);
  if(!type)
    return -1;
/*
 * Establish the latest clean model so that we have an up-to-date
 * representation of the UV model.
 */
  if(keep_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Determine statistics.
 */
  if(moddif(vlbob, &md, invpar.uvmin, invpar.uvmax))
    return -1;
/*
 * Assign the requested type as the function return value.
 */
  switch(type->id) {
  case UVRMS:
    *FLTPTR(outvals) = md.rms;
    break;
  case UVCHI:
    *FLTPTR(outvals) = md.chisq / md.ndata;
    break;
  case UVNVIS:
    *FLTPTR(outvals) = (md.ndata / 2L);
    break;
  };
  return no_error;
}

/*.......................................................................
 * Return image plane stats.
 *
 * Input:
 *  type     char *  The type of statistic to be returned.
 *                    rms    -  Map RMS flux.
 *                    mean   -  Map mean flux.
 *                    flux   -  Map total flux.
 *                    noise  -  Estimated noise.
 *                    bmin   -  Estimated minor axis of beam.
 *                    bmaj   -  Estimated major axis of beam.
 *                    bpa    -  Estimated major-axis position angle of beam.
 * Output:
 *  return  float    The requested statistic.
 */
static Template(imstat_fn)
{
/*
 * List names and enumeration identifiers for handled statistics.
 */
  enum {IMRMS, IMMEAN, IMNOISE, IMBMIN, IMBMAJ, IMBPA, IMDX, IMDY,
        IMDU, IMDV, IMNX, IMNY};
  static Enumpar types[] = {
    {"rms",IMRMS}, {"mean",IMMEAN}, {"noise",IMNOISE},
    {"bmin",IMBMIN}, {"bmaj",IMBMAJ}, {"bpa",IMBPA},
    {"dx", IMDX}, {"dy", IMDY}, {"du", IMDU}, {"dv", IMDV},
    {"nx", IMNX}, {"ny", IMNY}
  };
  static Enumtab *typtab=NULL; /* Enumerator symbol table */
  char *typename;        /* Statistic type name */
  Enumpar *type;         /* Pointer to element of types[] */
/*
 * Create the enumeration symbol table if necessary.
 */
  if(!typtab &&
     !(typtab=new_Enumtab(types, sizeof(types)/sizeof(Enumpar),"Image statistic")))
    return -1;
/*
 * Do we have a map to take parameters from?
 */
  if(nomap("imstat"))
    return -1;
/*
 * Invert the map if it is not ready for use?
 */
  if(vlbmap->domap && vlbmap->domap != MAP_IS_CLEAN &&
     invert_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Get the user's choice of statistic.
 */
  typename = *STRPTR(invals[0]);
/*
 * Determine the type of statistic member required.
 */
  type = find_enum(typtab, typename);
  if(!type)
    return -1;
/*
 * Assign the requested type as the function return value.
 */
  switch(type->id) {
  case IMRMS:
    *FLTPTR(outvals) = vlbmap->maprms;
    break;
  case IMMEAN:
    *FLTPTR(outvals) = vlbmap->mapmean;
    break;
  case IMNOISE:
    *FLTPTR(outvals) = vlbmap->noise;
    break;
  case IMBMIN:
    *FLTPTR(outvals) = radtoxy(vlbmap->e_bmin);
    break;
  case IMBMAJ:
    *FLTPTR(outvals) = radtoxy(vlbmap->e_bmaj);
    break;
  case IMBPA:
    *FLTPTR(outvals) = vlbmap->e_bpa * rtod;
    break;
  case IMDX:
    *FLTPTR(outvals) = radtoxy(vlbmap->xinc);
    break;
  case IMDY:
    *FLTPTR(outvals) = radtoxy(vlbmap->yinc);
    break;
  case IMDU:
    *FLTPTR(outvals) = vlbmap->uinc;
    break;
  case IMDV:
    *FLTPTR(outvals) = vlbmap->vinc;
    break;
  case IMNX:
    *FLTPTR(outvals) = (float) vlbmap->nx;
    break;
  case IMNY:
    *FLTPTR(outvals) = (float) vlbmap->ny;
    break;
  };
  return no_error;
}

/*.......................................................................
 * Plot amplitude versus projected UV radius.
 *
 * Input:
 *  pa     float   The UV plane position angle to project the data onto.
 *  telnam  char * The name of a telescope to highlight.
 *  uvmin  float   The min UV radius to display (user UV distance units).
 *  uvmax  float   The max UV radius (or 0 for full range) to plot.
 *                 (User UV distance units).
 *  ampmin float   The min amplitude to be displayed.
 *  ampmax float   The max amplitude to be displayed.
 *  phsmin float   The min phase to be displayed (degrees).
 *  phsmax float   The max phase to be displayed (degrees).
 *  docur logical  Default TRUE. If false, bypass interactive mode.
 */
static Template(uvprj_fn)
{
  Telspec *ts = NULL;/* Telescope specfication descriptor */
  char *opts="m3";   /* Default plot options */
  int docur=1;       /* Default to use of cursor where available */
  float uvmin=0.0f;  /* Default to full UV range */
  float uvmax=0.0f;
  float ampmin=0.0f; /* Default to full amplitude range */
  float ampmax=0.0f;
  float phsmin=0.0f; /* Default to full phase range */
  float phsmax=0.0f;
  float pa=0.0f;     /* Position angle (radians) */
  int modified;      /* Returned true by uvradplt() if data were modified */
  int iret;          /* Return value from uvradplt() */
/*
 * Make sure that we have data to be plotted.
 */
  if(nodata("projplot", OB_SELECT))
    return -1;
/*
 * Selectively overide defaults where given. NB. Fallthrough is deliberate.
 */
  switch(npar) {
  case 9:
    docur = *LOGPTR(invals[8]);
  case 8:
    phsmax = *FLTPTR(invals[7]) * dtor;
  case 7:
    phsmin = *FLTPTR(invals[6]) * dtor;
  case 6:
    ampmax = *FLTPTR(invals[5]);
  case 5:
    ampmin = *FLTPTR(invals[4]);
  case 4:
    uvmax = uvtowav(*FLTPTR(invals[3]));
  case 3:
    uvmin = uvtowav(*FLTPTR(invals[2]));
  case 2:
    ts = read_Telspec(vlbob, *STRPTR(invals[1]), NULL, 0);
    if(!ts)
      return -1;
  case 1:
    pa = *FLTPTR(invals[0]) * dtor;
  };
/*
 * Establish the latest clean model in order to be able to display it.
 */
  if(keep_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Open the display device.
 */
  if(make_open() == -1)
    return -1;
/*
 * Does the pflags options string contain any entries?
 * If so override the default options string.
 */
  if(**STRPTR(&pflags)=='\0')
    lprintf(stdout, "Using default options string \"%s\"\n", opts);
  else {
    opts = *STRPTR(&pflags);
    lprintf(stdout,
       "Overriding default options with user defined pflags=\"%s\"\n", opts);
  };
/*
 * Plot the data.
 */
  iret = uvradplt(vlbob, ts, docur, opts, 1, pa, uvmin, uvmax,
		  ampmin, ampmax, phsmin, phsmax, &modified);
/*
 * If the data were modified, mark the map and beam as out of date.
 */
  if(vlbmap && modified)
    vlbmap->domap = vlbmap->dobeam = MAP_IS_STALE;
  return iret ? -1 : no_error;
}

/*.......................................................................
 * Append normal models to continuum model.
 */
static Template(setcont_fn)
{
  int ncmp;  /* The number of components to be transfered */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("setcont", OB_INDEX))
    return -1;
/*
 * Count the number of components that are to be appended.
 */
  ncmp = vlbob->model->ncmp + vlbob->newmod->ncmp;
/*
 * Make sure that there are components to appended, so that the
 * map is not marked as invalid unnecessarily.
 */
  if(ncmp>0) {
/*
 * The UV data will have to be re-inverted to produce the modified
 * residual map. (No need to flag the beam for re-calculation).
 */
    if(vlbmap)
      vlbmap->domap = MAP_IS_STALE;
/*
 * Append the tentative and established models to the continuum models.
 */
    if(setcmod(vlbob, 1))
      return -1;
  };
/*
 * Inform user.
 */
  lprintf(stdout, "Added %d components to the continuum model.\n", ncmp);
  return no_error;
}

/*.......................................................................
 * Set the color-map parameters for the next maplot.
 *
 * Input:
 *  type     char *   The type of color-map none|grey|color|etc....
 *  contra   float  Contrast of color ramp (normally 1.0).
 *  bright   float  Brightness at the center color index (normally 0.5).
 */
static Template(mapcol_fn)
{
  Ctable *ctab = mappar.ctab;
/*
 * Install any new selections.
 */
  switch(npar) {
  default:
  case 3:
    ctab->bright = *FLTPTR(invals[2]);
  case 2:
    ctab->contra = *FLTPTR(invals[1]);
  case 1:
    if(!get_Cmap(ctab, *STRPTR(invals[0])))
      return -1;
  case 0:
    break;
  };
/*
 * Report the current settings.
 */
  lprintf(stdout, "Mapplot colormap: %s, contrast: %g brightness: %g.\n",
	  ctab->cmap->name, ctab->contra, ctab->bright);
  return no_error;
}

/*.......................................................................
 * Set the color-map transfer-function parameters for the next maplot.
 *
 * Input:
 *  type     char *   The type of color-map linear|logarithmic.
 *  vmin    float     Min data value to display.
 *  vmax    float     Max data value to display.
 */
static Template(mapfun_fn)
{
  Ctable *ctab = mappar.ctab;
/*
 * Install any new selections.
 */
  switch(npar) {
  default:
  case 3:
    ctab->vmax = *FLTPTR(invals[2]);
  case 2:
    ctab->vmin = *FLTPTR(invals[1]);
  case 1:
    ctab->tran = get_Cmtran(*STRPTR(invals[0]));
  case 0:
    break;
  };
/*
 * Report the current settings.
 */
  lprintf(stdout, "Mapplot transfer-function = %s, Data range = ",
	  name_Cmtran(ctab->tran));
  if(fabs(ctab->vmax - ctab->vmin) < 1.0e-15)
    lprintf(stdout, "data min -> data max.\n");
  else
    lprintf(stdout, "%g -> %g Jy.\n", ctab->vmin, ctab->vmax);
  return no_error;
}

/*.......................................................................
 * Set the mapplot clean beam placement.
 *
 * Input:
 *  xcenter    float   The fractional X-axis position of the beam center.
 *  ycenter    float   The fractional Y-axis position of the beam center.
 *  minsize    float   The min beam size wrt the size of the plot.
 *  maxsize    float   The max beam size wrt the size of the plot.
 */
static Template(beamloc_fn)
{
  MaplotBeam *mpb = &mappar.mpb;
/*
 * Install any new selections.
 */
  switch(npar) {
  default:
  case 4:
    mpb->maxsize = fabs(*FLTPTR(invals[3]));
  case 3:
    mpb->minsize = fabs(*FLTPTR(invals[2]));
  case 2:
    mpb->yc = *FLTPTR(invals[1]);
  case 1:
    mpb->xc = *FLTPTR(invals[0]);
  case 0:
    break;
  };
/*
 * Report the current settings.
 */
  if(mpb->xc < 0.0f || mpb->xc > 1.0f ||
     mpb->yc < 0.0f || mpb->yc > 1.0f) {
    lprintf(stdout, "Mapplot will not plot a clean beam ellipse.\n");
  } else {
    lprintf(stdout, 
           "Mapplot clean beam ellipse center: %g,%g. Size range: %g -> %g.\n",
	    mpb->xc, mpb->yc, mpb->minsize, mpb->maxsize);
  };
  return no_error;
}

/*.......................................................................
 * Set the mapplot polarization vector display parameters.
 *
 * Input:
 * scale float   The length in radians of a polarization vector of 1.0
 *               map intensity units.
 * icut  float   The minimum value of the normal intensity map at which
 *               polarization vectors should be drawn.
 * pcut  float   The minimum value of the polarization intensity map at
 *               which polarization vectors should be drawn.
 * dx,dy   int   Vectors will be drawn at every dx'th X-axis cell and
 *               every dy'th Y-axis cell.
 */
static Template(polvec_fn)
{
  MaplotVect *vect = &mappar.vect;
/*
 * Install any new selections.
 */
  switch(npar) {
  default:
  case 5:
    vect->dy = *INTPTR(invals[4]);
  case 4:
    vect->dx = *INTPTR(invals[3]);
  case 3:
    vect->pcut = *FLTPTR(invals[2]);
  case 2:
    vect->icut = *FLTPTR(invals[1]);
  case 1:
    vect->scale = xytorad(fabs(*FLTPTR(invals[0])));
  case 0:
    break;
  };
/*
 * Limit values.
 */
  if(vect->dx < 1)
    vect->dx = 1;
  if(vect->dy < 1)
    vect->dy = 1;
/*
 * Report the current settings.
 */
  lprintf(stdout, "Give polarization vectors lengths of %g %s/Jy.\n",
	  radtoxy(vect->scale), mapunits(U_TLAB));
  lprintf(stdout,
    "Draw vectors where unpolarized flux > %g Jy and polarized flux > %g Jy.\n",
	  vect->icut, vect->pcut);
  lprintf(stdout, "Draw polarization vectors in every ");
  if(vect->dx != 1)
    lprintf(stdout, "%d%s ", vect->dx, ordinal_suffix(vect->dx));
  lprintf(stdout, "X pixel and in every ");
  if(vect->dy != 1)
    lprintf(stdout, "%d%s ", vect->dy, ordinal_suffix(vect->dy));
  lprintf(stdout, "Y pixel.\n");
  return no_error;
}

/*.......................................................................
 * Display the state of the configuration parameters by writing the
 * parameter command file that would be written by save, but to the
 * screen.
 */
static Template(showpar_fn)
{
/*
 * The parameter setup is meaningless without an observation.
 */
  if(nodata("showpar", OB_INDEX))
    return -1;
  lprintf(stdout, "Difmap configuration state:\n");
  return wrtpars(NULL, NULL);
}

/*.......................................................................
 * Plot visibility spectra from the current observation.
 *
 * Input:
 *  avmode      char * The visibility averaging mode.
 *  ca,cb        int   Channel range.
 *  amin,amax  float   Amplitude plot range.
 *  pmin,pmax  float   Phase plot range.
 *  npage        int   The max number of pages to be plotted when plotting
 *                     to a non-interactive device.
 */
static Template(specpl_fn)
{
  Enumpar *epar;   /* A symbol table entry */
  int ca, cb;      /* Channel range */
  float amin,amax; /* Amplitude range */
  float pmin,pmax; /* Phase range */
  SpAvMode avmode; /* The visibility time-averaging mode */  
  int npage = 0;   /* Non-interactive page limit if non-zero */
  int nplot;
/*
 * Make sure that an observation has been read.
 */
  if(nodata("specplot", OB_INDEX))
    return -1;
/*
 * Make sure that a PGPLOT device is open.
 */
  if(make_open() == -1)
    return -1;
/*
 * Set defaults.
 */
  ca = cb = -1;
  amin = amax = 0.0f;
  pmin = pmax = 0.0f;
  avmode = SP_VECTOR;
  nplot = 0;
/*
 * Get the arguments.
 */
  switch(npar) {
  case 9:
    npage = *INTPTR(invals[8]);
  case 8:
    pmax = *FLTPTR(invals[7]);
  case 7:
    pmin = *FLTPTR(invals[6]);
  case 6:
    amax = *FLTPTR(invals[5]);
  case 5:
    amin = *FLTPTR(invals[4]);
  case 4:
    cb = *INTPTR(invals[3]) - 1;
  case 3:
    ca = *INTPTR(invals[2]) - 1;
  case 2:
    if((epar=find_enum(vlbspec->avsym, *STRPTR(invals[1]))) == NULL)
      return -1;
    avmode = (SpAvMode) epar->id;
  case 1:
    nplot = *INTPTR(invals[0]);
  };
/*
 * Install axis ranges.
 */
  if(sp_set_axes(vlbob, vlbspec, ca, cb, amin, amax, pmin, pmax))
    return -1;
/*
 * Initialize the averaging mode.
 */
  if(sp_set_options(vlbspec, nplot, vlbspec->xunit, avmode))
    return -1;
/*
 * Plot the selected spectra.
 */
  if(specplot(vlbob, vlbspec, 1, npage))
    return -1;
  return no_error;
}

/*.......................................................................
 * Set the list of baseline groups for a future invokation of specplot().
 *
 * Input:
 *  mode    char *   Baseline selection mode, from:
 *                    "split"   -  Show each baseline of the first group.
 *                    "group"   -  Show each baseline group as a whole.
 *  bstr    char *   One baseline selection list as a string.
 *  ...
 */
static Template(specb_fn)
{
  int i;
/*
 * Make sure that an observation has been read.
 */
  if(nodata("specbase", OB_INDEX))
    return -1;
/*
 * Parse arguments.
 */
  if(npar > 0) {
    Bgrplist *bgl = NULL; /* Baseline group. */
    SpBMode bmode;        /* Baseline selection mode */
    Enumpar *epar;        /* An enumeration symbol-table entry */
/*
 * Look up the baseline selection mode.
 */
    if((epar=find_enum(vlbspec->bmsym, *STRPTR(invals[0]))) == NULL)
      return -1;
    bmode = (SpBMode) epar->id;
/*
 * Get the list of baseline selection lists.
 */
    if(npar > 1) {
/*
 * Create an empty list of baseline groups.
 */
      bgl = new_Bgrplist();
      if(!bgl)
	return -1;
/*
 * Parse each baseline group argument and append it to the bgl list.
 */
      for(i=1; i<npar; i++) {
	if(!add_Basegrp(vlbob, bgl, NULL, *STRPTR(invals[i]))) {
	  bgl = del_Bgrplist(bgl);
	  return -1;
	};
      };
    };
/*
 * If new parameters were specified, record them in vlbspec.
 */
    if(sp_set_bgl(vlbob, vlbspec, bmode, bgl))
      return -1;
  };
/*
 * Display the contents of the current list.
 */
  {
    char awrk[80];
    Basegrp *bgrp;
    lprintf(stdout, "Specplot will plot baseline%s",
	    vlbspec->bmode==SP_GROUP ? " groups:\n" : "s of");
    for(bgrp=vlbspec->bgl->bgrp; bgrp; bgrp=bgrp->next) {
      if(write_Basegrp(vlbob, bgrp, sizeof(awrk), awrk) < 1)
	strcpy(awrk, "(Specification too long to display)");
      lprintf(stdout, " %s\n", awrk);
    };
  };
  return no_error;
}

/*.......................................................................
 * Set the list of polarizations for a future invokation of specplot().
 *
 * Input:
 *  pstr    char *   One polarization name.
 *  ...
 */
static Template(specp_fn)
{
  int i;
/*
 * Make sure that an observation has been read.
 */
  if(nodata("specpol", OB_INDEX))
    return -1;
  if(npar) {
/*
 * Create an empty list of polarizations.
 */
    Pollist *pl = new_Pollist();
    if(!pl)
      return -1;
/*
 * Parse each polarization argument and append it to the pl list.
 */
    for(i=0; i<npar; i++) {
      if(!add_Polnode(vlbob, pl, Stokes_id(*STRPTR(invals[i])))) {
	pl = del_Pollist(pl);
	return -1;
      };
    };
/*
 * Record the new list in vlbspec.
 */
    if(sp_set_pol(vlbob, vlbspec, pl))
      return -1;
  };
/*
 * Display the contents of the current list.
 */
  lprintf(stdout, "Specplot polarization selections: ");
  if(vlbspec->pl) {
    Polnode *pn;
    for(pn=vlbspec->pl->head; pn; pn=pn->next)
      lprintf(stdout, "%s%s", Stokes_name(pn->pol), pn->next ? ", " : "\n");
  } else {
    lprintf(stdout, "(default)\n");
  };
  return no_error;
}

/*.......................................................................
 * Set the time range for a future invokation of specplot().
 *
 * Input:
 *  stime     char *  The start of the range as a ddd/hh:mm:ss string.
 *  etime     char *  The end of the range as a ddd/hh:mm:ss string.
 *  interval float    The repeat interval within the range.
 */
static Template(spect_fn)
{
  double stime;   /* Start time (seconds) */
  double etime;   /* End time (seconds) */
  double scan;    /* Scan length or interval (seconds) */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("spectime", OB_INDEX))
    return -1;
/*
 * Get the start time.
 */
  if(npar > 0) {
/*
 * Get the start time.
 */
    if((stime=read_ut(*STRPTR(invals[0]), NULL)) < 0)
      return -1;
/*
 * Get the end time.
 */
    if(npar > 1) {
      if((etime=read_ut(*STRPTR(invals[1]), NULL)) < 0)
	return -1;
    } else {
      etime = vlbob->rec[vlbob->nrec-1].integ->ut;
    };
/*
 * Read the scan time.
 */
    if(npar > 2) {
      scan = *FLTPTR(invals[2]) * 60.0; /* Convert from minutes to seconds */
    } else {
      scan = fabs(etime - stime);
    };
/*
 * Install the new time range.
 */
    if(sp_set_times(vlbob, vlbspec, stime, etime, scan))
      return -1;
  };
/*
 * Report the current time range setting.
 */
  {
    char awrk[80];
    write_ut(vlbspec->stime, sizeof(awrk), awrk);
    lprintf(stdout, "Specplot time range %s - ", awrk);
    write_ut(vlbspec->etime, sizeof(awrk), awrk);
    lprintf(stdout, "%s,  scan %s=%g mins\n", awrk,
	    vlbspec->scan<0.0 ? "separation" : "length",
	    fabs(vlbspec->scan/60.0f));
  };
  return no_error;
}

/*.......................................................................
 * Set the UV range attributes for a future invokation of specplot().
 *
 * Input:
 *  uvmin     float   The min UV radius (wavelengths).
 *  uvmax     float   The max UV radius, or 0 for the max available.
 *  uvstep    float   The step size to break uvmin -> uvmax into, or
 *                    0 for uvmax - uvmin.
 */
static Template(specuv_fn)
{
  float uvmin;  /* Min UV radius */
  float uvmax;  /* Max UV radius */
  float uvstep; /* Iterator step size */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("spectime", OB_INDEX))
    return -1;
/*
 * Install new values.
 */
  if(npar > 0) {
    uvmin = uvmax = uvstep = 0.0f;
    switch(npar) {
    case 3:
      uvstep = uvtowav(*FLTPTR(invals[2]));
    case 2:
      uvmax = uvtowav(*FLTPTR(invals[1]));
    case 1:
      uvmin = uvtowav(*FLTPTR(invals[0]));
    };
    if(sp_set_uvrange(vlbob, vlbspec, uvmin, uvmax, uvstep))
      return -1;
  };
/*
 * Report the current UV range.
 */
  uvmin = wavtouv(vlbspec->uvr.uvmin);
  uvmax = wavtouv(vlbspec->uvr.uvmax);
  uvstep = wavtouv(vlbspec->uvr.uvstep);
  lprintf(stdout, "Specplot UV range: uvmin=%g  uvmax=%g  uvstep=%g (%s)\n",
	  uvmin, uvmax, uvstep, uvwunits(U_TLAB));
  return no_error;
}

/*.......................................................................
 * Set the smoothing options for a future invokation of specplot().
 *
 * Input:
 *  xtype      char *   The units of the fwhm of the smoothing function.
 *                       channels, frequency
 *  stype      char *   The type of smoothing required from:
 *                       raw, hanning, gaussian, boxcar, sinc.
 *  fwhm      float     The FWHM of the smoothing function.
 */
static Template(specsm_fn)
{
  Enumpar *epar;  /* A symbol table entry */
  SpXunit xunit;  /* The enumeration identifier of the chosen axis type */
  SmType smtype;  /* The enumeration identifier of the chosen smoothing type */
  float fwhm;     /* Full-width-half-maximum of smoothing function */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("specsmooth", OB_INDEX))
    return -1;
/*
 * Set defaults.
 */
  xunit = SP_CHAN;
  smtype = SM_NONE;
  fwhm = 0.0f;
/*
 * Get the X-axis type.
 */
  switch(npar) {
  case 3:
    fwhm = *FLTPTR(invals[2]);
    if(fwhm < 0.0f) {
      lprintf(stderr, "Unacceptable fwhm: %g\n", fwhm);
	return -1;
    };
  case 2:
    if((epar=find_enum(vlbspec->smsym, *STRPTR(invals[1]))) == NULL)
      return -1;
    smtype = (SmType) epar->id;
  case 1:
    if((epar=find_enum(vlbspec->xtsym, *STRPTR(invals[0]))) == NULL)
      return -1;
    xunit = (SpXunit) epar->id;
  };
/*
 * Record the new defaults.
 */
  if(npar > 0 && sp_set_smooth(vlbspec, xunit, smtype, fwhm))
    return -1;
/*
 * Report the current settings.
 */
  lprintf(stdout, "Specplot smoothing:  units=%s  window=%s  fwhm=%g\n",
	  name_enum(vlbspec->xtsym, vlbspec->xunit, "(unknown)"),
	  name_enum(vlbspec->smsym, vlbspec->smooth.type, "(unknown)"),
	  vlbspec->smooth.fwhm);
  return no_error;
}

/*.......................................................................
 * Set specplot plot modes.
 *
 * Input:
 *  nplot        int    The number of plots per page (0 -> 3).
 *  xunit       char *  The units along the X axis.
 *                       channels      - Channel indexes.
 *                       frequency     - Frequency.
 *  avmode      char *  Visibility time-averaging mode.
 *                       vector        - Vector average.
 *                       scalar        - Scalar average.
 *  flags       char *  An array of option flags.
 */
static Template(specop_fn)
{
  Enumpar *epar;   /* A symbol table entry */
  SpXunit xunit;   /* The enumeration identifier of the chosen axis type */
  char *flags;     /* Array of flags */
/*
 * Make sure that an observation has been read.
 */
  if(nodata("specopt", OB_INDEX))
    return -1;
/*
 * Set defaults.
 */
  xunit = vlbspec->xunit;
  flags = NULL;
/*
 * Get the command-line arguments.
 */
  switch(npar) {
  case 2:
    flags = *STRPTR(invals[1]);
  case 1:
    if((epar=find_enum(vlbspec->xtsym, *STRPTR(invals[0]))) == NULL)
      return -1;
    xunit = (SpXunit) epar->id;
  default:
  case 0:
    break;
  };
/*
 * Record the new modes.
 */
  if(npar > 0 && sp_set_options(vlbspec, vlbspec->nplot, xunit, SP_VECTOR))
    return -1;
/*
 * Install and/or get the current set of option flags.
 */
  flags = sp_set_flags(vlbspec, flags);
  if(!flags)
    return -1;
/*
 * Report the current settings.
 */
  lprintf(stdout, "Specplot x-axis=%s  flags=\"%s\".\n",
	  name_enum(vlbspec->xtsym, vlbspec->xunit, "(unknown)"), flags);
  return no_error;
}

/*.......................................................................
 * Set the order of specplot specification keys.
 *
 * Input:
 *  key1,...keyn  The names of specification keys in the desired sort
 *                order.
 */
static Template(specso_fn)
{
  Enumpar *epar;        /* A symbol table entry */
  SpKey keys[SP_NKEY];  /* The array of npar keys */
  int i;
/*
 * Make sure that an observation has been read.
 */
  if(nodata("specorder", OB_INDEX))
    return -1;
/*
 * Get the key names.
 */
  for(i=0; i<npar; i++) {
    if((epar=find_enum(vlbspec->keysym, *STRPTR(invals[i]))) == NULL)
      return -1;
    keys[i] = (SpKey) epar->id;
  };
/*
 * Record the new selection order.
 */
  if(npar > 0 && sp_set_order(vlbspec, keys, npar))
    return -1;
/*
 * Report the current settings.
 */
  lprintf(stdout, "Specplot selection order:");
  for(i=0; i<vlbspec->nkey; i++) {
    lprintf(stdout, " %s",
	    name_enum(vlbspec->keysym, vlbspec->key[i], "(unknown)"));
  };
  lprintf(stdout, "\n");
  return no_error;
}

/*.......................................................................
 * Return the value of the pixel nearest to a given center-relative
 * map coordinate.
 *
 * Input:
 *  x,y     float    The coordinate to return a value for.
 * Output:
 *  return  float    The corresponding pixel value.
 */
static Template(mapval_fn)
{
  float x,y;  /* The selected map coordinates */
  int ix,iy;  /* The index of the map pixel to probe */
/*
 * Do we have a map to take parameters from?
 */
  if(nomap("mapvalue"))
    return -1;
/*
 * Invert the map if it is not ready for use?
 */
  if(vlbmap->domap && vlbmap->domap != MAP_IS_CLEAN &&
     invert_fn(invals,0,outvals) == -1)
    return -1;
/*
 * Get the user's map coordinates.
 */
  x = *FLTPTR(invals[0]);
  y = *FLTPTR(invals[1]);
/*
 * Get the indexes of the pixel nearest to the given map coordinates.
 */
  ix = map_x_coord_to_pixel(vlbmap, xytorad(x));
  iy = map_y_coord_to_pixel(vlbmap, xytorad(y));
/*
 * Check bounds.
 */
  if(ix < vlbmap->maparea.ixmin || ix > vlbmap->maparea.ixmax ||
     iy < vlbmap->maparea.iymin || iy > vlbmap->maparea.iymax) {
    lprintf(stderr, "mapvalue: Coordinates out of range.\n");
    return -1;
  };
/*
 * Return the map pixel value.
 */
  *FLTPTR(outvals) = vlbmap->map[ix + vlbmap->nx * iy];
  return no_error;
}

/*.......................................................................
 * Shift the phase center of the observation to a given RA,DEC position.
 *
 * Input:
 *  ra    char *   The Right Ascension expressed in sexagesimal hours.
 *  dec   char *   The Declination expressed in sexagesimal degrees.
 */
static Template(shiftto_fn)
{
  char *ra_s,*dec_s; /* The strings that contain the target RA,Dec */
  double ra,dec;     /* The target Right Ascension and Declination (radians) */
  double north,east; /* The direction-cosine northward and eastward shifts */
/*
 * The following objects are used to pass arguments to shift_fn().
 */
  static float n_shift,e_shift;     
  static Descriptor d_north  = {'f',0,R_ONLY,1,{1,1,1},&n_shift};
  static Descriptor d_east  = {'f',0,R_ONLY,1,{1,1,1},&e_shift};
  static Descriptor *shift_args[3] = {&d_east, &d_north};
/*
 * Check that there is something to shift.
 */
  if(nodata("shiftto", OB_INDEX))
    return -1;
/*
 * Get the arguments that contain the RA and Dec as sexagesimal strings.
 */
  ra_s = *STRPTR(invals[0]);
  dec_s = *STRPTR(invals[1]);
/*
 * Parse the Right Ascension argument.
 */
  if(parse_sexagesimal_string(ra_s, &ra, NULL)) {
    lprintf(stderr, "shiftto: Bad Right Ascension: %s\n", ra_s);
    return -1;
  };
  ra *= htor;
/*
 * Parse the Declination argument.
 */
  if(parse_sexagesimal_string(dec_s, &dec, NULL)) {
    lprintf(stderr, "shiftto: Bad Declination: %s\n", dec_s);
    return -1;
  };
  dec *= dtor;
/*
 * Compute the direction cosine shifts required to shift from the
 * recorded phase center to the specified position, taking account
 * of any shifts that have already been applied.
 */
  east = radec_to_l(vlbob->source.ra, vlbob->source.dec, ra, dec, vlbob->proj)
    - vlbob->geom.east;
  north = radec_to_m(vlbob->source.ra, vlbob->source.dec, ra, dec, vlbob->proj)
    - vlbob->geom.north;
/*
 * Let the shift command perform the actual shift.
 */
  n_shift = radtoxy(-north);
  e_shift = radtoxy(-east);
  return shift_fn(shift_args, 2, outvals);
}

/*.......................................................................
 * Create polarization intensity and angle maps above and below the
 * map of the current selection. Note that the beam array is used as
 * a work area during this operation, so vlbmap->dobeam will be 1 after
 * this function has been called.
 *
 * Input:
 *  docln    int    If true attempt to restore each map.
 * Output:
 *  return   int    0 - OK.
 *                 -1 - Error.
 */
static int make_polmap(int docln)
{
  int xa,xb;     /* The delimiting X-axis indexes of the inner quarter of */
                 /*  the map. */
  int ya,yb;     /* The delimiting Y-axis indexes of the inner quarter of */
                 /*  the map. */
  int ix,iy;     /* Map x,y indexes */
  float *mptr;   /* A pointer into vlbmap->map[] */
  float *uptr;   /* A pointer into the U map */
  float *qptr;   /* A pointer into the Q map */
  float *magptr; /* A pointer into the polarization magnitude map */
  float *angptr; /* A pointer into the polarization angle map */
  Stokes pol;    /* The currently selected polarization */
  int i;
/*
 * Check that a map and beam have been allocated and that an observation
 * has been read.
 */
  if(nomap("make_polmap") || nodata("make_polmap", OB_SELECT))
    return -1;
/*
 * If making a restored map we need different models for different
 * polarizations.
 */
  if(docln && !multi_model_mode) {
    lprintf(stderr,
         "Models are needed for each polarization - see help multi_model.\n");
    return -1;
  };
/*
 * Don't remake the map if it is already in the map array.
 */
  if(vlbmap->domap == (docln ? MAP_IS_PCLN : MAP_IS_PMAP))
    return no_error;
/*
 * We will need both U and Q in order to make the polarization
 * map.
 */
  if(get_Obpol(vlbob, SU, 0, NULL) || get_Obpol(vlbob, SQ, 0, NULL)) {
    lprintf(stderr,
	 "To make polarization maps you need either U and Q, or LR and RL.\n");
    return -1;
  };
/*
 * Keep a record of the currently selected polarization.
 */
  pol = vlbob->stream.pol.type;
/*
 * If displaying a restored map, make sure that models exist for all of
 * the specified polarizations.
 */
  if(docln) {
    Stokes bad_pol[3];   /* The 'nbad' polarizations which don't have models */
    int nbad = 0;
/*
 * Does the currently selected polarization have a model?
 */
    if(vlbob->model->ncmp + vlbob->newmod->ncmp < 1)
      bad_pol[nbad++] = pol;
/*
 * Do we have a model for Q?
 */
    if(pol!=SQ && !have_ModelEntry(vlbob->mtab, vlbob->stream.cl, SQ, 1))
      bad_pol[nbad++] = SQ;
/*
 * Do we have a model for U?
 */
    if(pol!=SU && !have_ModelEntry(vlbob->mtab, vlbob->stream.cl, SU, 1))
      bad_pol[nbad++] = SU;
/*
 * Report the lack of any of the models.
 */
    if(nbad) {
      lprintf(stderr, "Please clean or modelfit polarization%s:",
	      nbad==1 ? "" : "s");
      for(i=0; i<nbad; i++)
	lprintf(stderr, " %s", Stokes_name(bad_pol[i]));
      lprintf(stderr, ".\n");
      return -1;
    };
  };
/*
 * Determine the start and end indexes of the central map.
 */
  xa = vlbmap->nx/4;
  ya = vlbmap->ny/4;
  xb = 3*xa - 1;
  yb = 3*ya - 1;
/*
 * If making a restored map we need to estimate the size of the restoring
 * beam from the dirty beam.
 */
  if(docln && vlbmap->dobeam && invert_fn(NULL, 0, NULL))
    return -1;
/*
 * We will be using the beam array as working storage, so we don't
 * want 'uvinvert' to overwrite it with an update of the beam.
 */
  vlbmap->dobeam = 0;
/*
 * Make a Q map.
 */
  if(pol != SQ) {
    if(ob_select(vlbob, 0, vlbob->stream.cl, SQ))
      return -1;
    vlbmap->domap = MAP_IS_STALE;
  };
/*
 * Create the Q map.
 */
  if(vlbmap->domap && !(docln && vlbmap->domap==MAP_IS_CLEAN) &&
     invert_fn(NULL, 0, NULL))
    return polmap_error(pol);
/*
 * Restore the Q map if needed
 */
  if(docln && (!vlbmap->ncmp || vlbmap->domap!=MAP_IS_CLEAN) &&
     restore_fn(NULL,0,NULL) == -1)
    return polmap_error(pol);
/*
 * Copy the significant inner quarter of the Q map into the bottom quarter
 * of the beam array.
 */
  qptr = vlbmap->beam;
  mptr = vlbmap->map + xa + ya * vlbmap->nx;
  for(iy=ya; iy<=yb; iy++,mptr+=vlbmap->nx/2) {
    for(ix=xa; ix<=xb; ix++)
      *qptr++ = *mptr++;
  };
/*
 * Get the U stokes parameter data.
 */
  if(ob_select(vlbob, 0, vlbob->stream.cl, SU))
    return polmap_error(pol);
/*
 * Create the residual U map.
 */
  vlbmap->domap = MAP_IS_STALE;
  if(invert_fn(NULL, 0, NULL))
    return polmap_error(pol);
/*
 * Create the restored U map if needed.
 */
  if(docln && restore_fn(NULL,0,NULL) == -1)
    return polmap_error(pol);
/*
 * Copy the significant inner quarter of the U map into the second quarter
 * of the beam array.
 */
  uptr = vlbmap->beam + vlbmap->ny/4 * vlbmap->nx;
  mptr = vlbmap->map + xa + ya * vlbmap->nx;
  for(iy=ya; iy<=yb; iy++,mptr+=vlbmap->nx/2) {
    for(ix=xa; ix<=xb; ix++)
      *uptr++ = *mptr++;
  };
/*
 * Create a map of the stokes parameter that was selected on entry to
 * this command.
 */
  if(pol != vlbob->stream.pol.type) {
    if(ob_select(vlbob, 0, vlbob->stream.cl, pol))
      return polmap_error(pol);
/*
 * Create the residual map.
 */
    vlbmap->domap = MAP_IS_STALE;
    if(invert_fn(NULL, 0, NULL))
      return polmap_error(pol);
/*
 * Create the restored map if possible.
 */
    if(docln && restore_fn(NULL,0,NULL) == -1)
      return polmap_error(pol);
  };
/*
 * Convert the Q and U maps into polarized intensity and polarized
 * angle images, placing them in the margins just above and below
 * the image that is currently in the map array.
 */
  qptr = vlbmap->beam;
  uptr = vlbmap->beam + vlbmap->ny/4 * vlbmap->nx;
  magptr = vlbmap->map;
  angptr = vlbmap->map + 3*vlbmap->ny/4 * vlbmap->nx;
  for(iy=0; iy<vlbmap->ny/2; iy++) {
    for(ix=0; ix<vlbmap->nx/2; ix++) {
      float q = *qptr++;
      float u = *uptr++;
      *magptr++ = sqrt(q*q + u*u);
      *angptr++ = u==0.0 && q==0.0 ? 0.0 : 0.5 * atan2(u,q);
    };
  };
/*
 * The beam array will need to be recomputed before it can be used again.
 */
  vlbmap->dobeam = 1;
/*
 * Advertise the new contents of the map.
 */
  vlbmap->domap = docln ? MAP_IS_PCLN : MAP_IS_PMAP;
  return no_error;
}

/*.......................................................................
 * This is the private error return function of make_polmap().
 *
 * Input:
 *  pol   Stokes  The polarization that was selected when make_polmap()
 *                was invoked.
 * Output:
 *  return   int  The error-return code of make_polmap().
 */
static int polmap_error(Stokes pol)
{
/*
 * Mark both the map and the beam as needing to be updated.
 */
  vlbmap->domap = vlbmap->dobeam = MAP_IS_STALE;
/*
 * Attempt to select the polarization that was selected when make_polmap()
 * was first entered.
 */
  (void) ob_select(vlbob, !multi_model_mode, vlbob->stream.cl, pol);
/*
 * Return the error code of make_polmap().
 */
  return -1;
}

/*.......................................................................
 * Write the models associated with all mapped sets of channel ranges
 * and polarizations to a file.
 *
 * Input:
 *  name   char *  The file name for the multi-model file.
 */
static Template(write_models_fn)
{
  char *filename;  /* The name to give the new file */
/*
 * Check that an observation has been read.
 */
  if(nodata("write_models", OB_INDEX))
    return -1;
/*
 * Get the filename.
 */
  filename = *STRPTR(invals[0]);
/*
 * Record the model of the current stream in the model table so that
 * it gets saved.
 */
  if(ob_ready(vlbob, OB_SELECT, NULL) && ob_record_select_model(vlbob))
    return -1;
/*
 * Write the models.
 */
  if(write_ModelTable(vlbob->mtab, filename))
    return -1;
  return no_error;
}

/*.......................................................................
 * Read models from a given multi-model file into the table of models
 * associated with different channel/polarization selections.
 *
 * Input:
 *  name   char *  The file name for the multi-model file.
 */
static Template(read_models_fn)
{
  char *filename;  /* The name of the input file */
/*
 * Check that an observation has been read.
 */
  if(nodata("read_models", OB_INDEX))
    return -1;
/*
 * Get the filename.
 */
  filename = *STRPTR(invals[0]);
/*
 * Write the models.
 */
  if(read_ModelTable(vlbob->mtab, filename))
    return -1;
/*
 * If we are in multi-model mode and a channel/polarization selection
 * has already been made, replace the current model with the appropriate
 * model from the modified table.
 */
  if(multi_model_mode && ob_ready(vlbob, OB_SELECT, NULL) &&
     ob_install_select_model(vlbob))
      return -1;
/*
 * If not in multi-model mode, warn the user that the models won't be
 * visible.
 */
  if(!multi_model_mode) {
    lprintf(stdout,
   "read_models: Warning: Multi-model mode disabled - see help multi_model.\n");
  };
  return no_error;
}

/*.......................................................................
 * Clear all models associated with multiple channel/polarization
 * selections.
 */
static Template(clear_models_fn)
{
/*
 * Check that an observation has been read.
 */
  if(nodata("clear_models", OB_INDEX))
    return -1;
/*
 * Clear the table of models.
 */
  if(clear_ModelTable(vlbob->mtab))
    return -1;
  return no_error;
}

/*.......................................................................
 * Enable or disable keeping different models for different selections.
 *
 * Input:
 *  domulti  logical  True to enable multi-model mode.
 */
static Template(multi_model_fn)
{
/*
 * Has the user requested a change in status?
 */
  if(npar > 0) {
    int was_multi = multi_model_mode;
    multi_model_mode = *LOGPTR(invals[0]);
/*
 * If we are switching into multi-model mode and a channel/polarization
 * selection has already been made, replace the current model with
 * the appropriate model from the modified table.
 */
    if(multi_model_mode) {
      if(!was_multi && vlbob && ob_ready(vlbob, OB_SELECT, NULL) &&
	 ob_install_select_model(vlbob))
	return -1;
      if(vlbmap)
	vlbmap->domap = MAP_IS_STALE;
/*
 * If we are switching out of multi-model mode, record the current model,
 * in the table of models.
 */
    } else {
      if(was_multi && vlbob && ob_ready(vlbob, OB_SELECT, NULL) &&
	 ob_record_select_model(vlbob))
	return -1;
    };
  };
/*
 * Report the current status.
 */
  if(multi_model_mode)
    lprintf(stdout,
	"Maintain separate models for each channel/polarization selection.\n");
  else
    lprintf(stdout,
	"Use one model for all channel/polarization selections.\n");
  return no_error;
}

/*.......................................................................
 * Add a marker to the list of markers that are to be drawn on subsequent
 * maps, specifying its position by its Right Ascension and Declination.
 *
 * Input:
 *  ra,dec     The Right Ascension and Declination (sexagesimal strings).
 *  symbol     The name of a marker symbol.
 *  color      The PGPLOT color index to use when displaying the marker.
 *  size       The character size to use relative to the normal size of 1.0.
 *  text       Optional anotation text to place next to the marker.
 *  just       The justification of the text (0 right, 0.5 center, 1 left).
 *  xpos       The x-offset of the justification point of the text from
 *             the marker (characters).
 *  ypos       The y-offset of the vertical center of the text from the
 *             the marker (characters).
 */
static Template(mark_radec_fn)
{
  float just = 0.0;     /* Default to right justification */
  float ypos = 0.0;     /* Default to having the text level with the marker */
  float xpos = 1.0;     /* Default to having the text start one character */
                        /*  width after the leftmost character. */
  int color = 11;       /* The color to use when plotting the marker */
  float size = 1.0;     /* The character size */
  char *text = NULL;    /* Optional annotation text */
  char *sym_s = "dot";  /* The name of the plot symbol */
  char *ra_s,*dec_s;    /* The RA and Dec as sexagesimal strings */
  double ra, dec;       /* The RA and Dec in radians */
  MarkerSymbol symbol;  /* The plot symbol code */
/*
 * Override the defaults with the user's arguments.
 */
  switch(npar) {
  case 9:
    ypos = *FLTPTR(invals[8]);
  case 8:
    xpos = *FLTPTR(invals[7]);
  case 7:
    just = *FLTPTR(invals[6]);
  case 6:
    text = *STRPTR(invals[5]);
  case 5:
    size = *FLTPTR(invals[4]);
  case 4:
    color = *INTPTR(invals[3]);
  case 3:
    sym_s = *STRPTR(invals[2]);
  default:
    dec_s = *STRPTR(invals[1]);
    ra_s = *STRPTR(invals[0]);
  };
/*
 * Parse the Right Ascension argument.
 */
  if(parse_sexagesimal_string(ra_s, &ra, NULL)) {
    lprintf(stderr, "mark_radec: Bad Right Ascension: %s\n", ra_s);
    return -1;
  };
  ra *= htor;
/*
 * Parse the Declination argument.
 */
  if(parse_sexagesimal_string(dec_s, &dec, NULL)) {
    lprintf(stderr, "mark_radec: Bad Declination: %s\n", dec_s);
    return -1;
  };
  dec *= dtor;
/*
 * Lookup the named symbol.
 */
  symbol = lookup_marker_symbol(mapmarkers, sym_s);
  if(symbol == MK_UNKNOWN)
    return -1;
/*
 * Add the marker definition to the list.
 */
  if(add_MarkerNode(mapmarkers, ra, dec, symbol, color, size, text, just,
		    xpos, ypos) == NULL)
    return -1;
  return no_error;
}

/*.......................................................................
 * Add a marker to the list of markers that are to be drawn on subsequent
 * maps, specifying its position as an X,Y position on the map.
 *
 * Input:
 *  x,y        The Right Ascension and Declination (sexagesimal strings).
 *  symbol     The name of a marker symbol.
 *  color      The PGPLOT color index to use when displaying the marker.
 *  size       The character size to use relative to the normal size of 1.0.
 *  text       Optional anotation text to place next to the marker.
 *  just       The justification of the text (0 right, 0.5 center, 1 left).
 *  xpos       The x-offset of the justification point of the text from
 *             the marker (characters).
 *  ypos       The y-offset of the vertical center of the text from the
 *             the marker (characters).
 */
static Template(mark_xy_fn)
{
  float just = 0.0;     /* Default to right justification */
  float ypos = 0.0;     /* Default to having the text level with the marker */
  float xpos = 1.0;     /* Default to having the text start one character */
                        /*  width after the leftmost character. */
  int color = 11;       /* The color to use when plotting the marker */
  float size = 1.0;     /* The character size */
  char *text = NULL;    /* Optional annotation text */
  char *sym_s = "dot";  /* The name of the plot symbol */
  float x,y;            /* The map X,Y coordinates of the source */
  double ra, dec;       /* The RA and Dec in radians */
  MarkerSymbol symbol;  /* The plot symbol code */
/*
 * We need a dataset in order to have a phase center to refer the shift
 * to.
 */
  if(nodata("mark_xy", OB_INDEX))
    return -1;
/*
 * Override the defaults with the user's arguments.
 */
  switch(npar) {
  case 9:
    ypos = *FLTPTR(invals[8]);
  case 8:
    xpos = *FLTPTR(invals[7]);
  case 7:
    just = *FLTPTR(invals[6]);
  case 6:
    text = *STRPTR(invals[5]);
  case 5:
    size = *FLTPTR(invals[4]);
  case 4:
    color = *INTPTR(invals[3]);
  case 3:
    sym_s = *STRPTR(invals[2]);
  default:
    y = xytorad(*FLTPTR(invals[1])) - vlbob->geom.north;
    x = xytorad(*FLTPTR(invals[0])) - vlbob->geom.east;
  };
/*
 * Convert from the specified x,y offset from the current center of the
 * map to Right Ascension and Declination.
 */
  ra = lmtora(vlbob->source.ra, vlbob->source.dec, x, y, vlbob->proj);
  dec = lmtodec(vlbob->source.ra, vlbob->source.dec, x, y, vlbob->proj);
/*
 * Lookup the named symbol.
 */
  symbol = lookup_marker_symbol(mapmarkers, sym_s);
  if(symbol == MK_UNKNOWN)
    return -1;
/*
 * Add the marker definition to the list.
 */
  if(add_MarkerNode(mapmarkers, ra, dec, symbol, color, size, text, just,
		    xpos, ypos) == NULL)
    return -1;
  return no_error;
}

/*.......................................................................
 * Clear the current list of markers.
 */
static Template(clear_markers_fn)
{
  clr_MarkerList(mapmarkers);
  lprintf(stdout, "Marker list cleared.\n");
  return no_error;
}

/*.......................................................................
 * Write the list of commands needed to reproduce the current list of
 * map markers, writing them either to a file or to stdout.
 *
 * Input:
 *  name   char *  The file name for the command file, or "" to
 *                 write to stdout.
 */
static Template(wmarkers_fn)
{
  FILE *fp = stdout;           /* Default to writing to stdout */
  char *fname = "stdout";      /* The name of the output file */
/*
 * Sanity check.
 */
  if(nodata("wmarkers", OB_INDEX))
    return -1;
/*
 * Get arguments.
 */
  switch(npar) {
  case 1:
/*
 * Override the stdout default?
 */
    if(*STRPTR(invals[0])[0] != '\0') {
/*
 * Get the model-file name.
 */
      fname = *STRPTR(invals[0]);
/*
 * Open the model file.
 */
      fp = fopen(fname, "w");
      if(fp == NULL) {
	lprintf(stderr, "wmarkers: Unable to create file: %s\n", fname);
	return -1;
      };
    };
  };
/*
 * Keep the user informed.
 */
  lprintf(stderr, "Writing marker commands to %s\n", fname);
/*
 * Compose the commands needed to restore each of the markers in
 * the list.
 */
  if(write_marker_commands(fp) || (fp!=stdout && fclose(fp)==EOF)) {
    lprintf(stderr, "wmarkers: Error writing file: %s\n", fname);
    return -1;
  };
  return no_error;
}

/*.......................................................................
 * Write the list of commands needed to restore the current list of
 * map markers, to a given stream.
 *
 * Input:
 *  fp    FILE *   The file to write to.
 * Output:
 *  return int     0 - OK.
 *                 1 - Error (unreported).
 */
static int write_marker_commands(FILE *fp)
{
  MarkerNode *marker;
  for(marker=mapmarkers->head; marker; marker=marker->next) {
    char ra_s[20];   /* The RA as a sexagesimal string */
    char dec_s[20];  /* The Dec as a sexagesimal string */
    if(lprintf(fp, "mark_radec %s, %s, %s, %d, %g, ",
	       sradhms(marker->ra, 3, 1, ra_s),
	       sraddms(marker->dec, 3, 1, dec_s),
	       lookup_marker_name(mapmarkers, marker->sym),
	       marker->color, marker->size) < 0)
      return 1;
    if(write_string_arg(fp, NULL, marker->text ? marker->text : ""))
      return 1;
    if(lprintf(fp, ", %g, %g, %g\n", marker->just, marker->xpos,
	       marker->ypos) < 0)
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Delete an Ra,Dec marker from the list of map markers.
 *
 * Input:
 *  ra,dec     The Right Ascension and Declination (sexagesimal strings).
 */
static Template(delmarker_fn)
{
  char *ra_s,*dec_s;    /* The RA and Dec as sexagesimal strings */
  double ra, dec;       /* The RA and Dec in radians */
  MarkerNode *nearest;  /* The marker to be deleted */
/*
 * Get the target Ra,Dec.
 */
  dec_s = *STRPTR(invals[1]);
  ra_s = *STRPTR(invals[0]);
/*
 * Parse the Right Ascension argument.
 */
  if(parse_sexagesimal_string(ra_s, &ra, NULL)) {
    lprintf(stderr, "delmarker: Bad Right Ascension: %s\n", ra_s);
    return -1;
  };
  ra *= htor;
/*
 * Parse the Declination argument.
 */
  if(parse_sexagesimal_string(dec_s, &dec, NULL)) {
    lprintf(stderr, "delmarker: Bad Declination: %s\n", dec_s);
    return -1;
  };
  dec *= dtor;
/*
 * Find the nearest marker node to the specified RA, Dec.
 */
  nearest = closest_MarkerNode(mapmarkers, ra, dec);
/*
 * No markers in list?
 */
  if(!nearest) {
    lprintf(stderr, "delmarker: No marker found.\n");
    return -1;
  };
/*
 * Delete the marker.
 */
  nearest = del_MarkerNode(mapmarkers, nearest);
  return no_error;
}

/*.......................................................................
 * Return an array of statistics for a given observable.
 *
 * Input:
 *  qty    char *   The name of an observable.
 * Output:
 *  return float[6] An array of statistics, containing:
 *                   nvis    - The number of visibilities used.
 *                   mean    - The mean of the observable.
 *                   error   - The standard error on the mean.
 *                   scatter - The scatter about the mean.
 *                   minval  - The minimum value of the observable.
 *                   maxval  - The maximum value of the observable.
 */
static Template(vis_stats_fn)
{
  VisStat stats;   /* The return container of ob_vis_stats() */
  float *fptr;     /* The return array */
  double cnvfac;   /* Conversion factor from internal to user units */
/*
 * List names for handled observables.
 */
  static Enumpar types[] = {
    {"amplitude", VS_AMP}, {"phase", VS_PHS},
    {"real", VS_REAL}, {"imaginary", VS_IMAG},
    {"umag", VS_UMAG}, {"vmag", VS_VMAG}, {"uvrad", VS_UVRAD}
  };
  static Enumtab *typtab=NULL; /* Enumerator symbol table */
  char *typename;              /* Observable type name */
  Enumpar *type;               /* Pointer to element of types[] */
/*
 * Create the enumeration symbol table if necessary.
 */
  if(!typtab &&
     !(typtab=new_Enumtab(types, sizeof(types)/sizeof(Enumpar),"Observable")))
    return -1;
/*
 * Make sure that we have data to examine.
 */
  if(nodata("vis_stats", OB_SELECT))
    return -1;
/*
 * Get the name of the observable.
 */
  typename = *STRPTR(invals[0]);
/*
 * Lookup the observable.
 */
  type = find_enum(typtab, typename);
  if(!type)
    return -1;
/*
 * Get the requested stats.
 */
  if(ob_vis_stats(vlbob, type->id, invpar.uvmin, invpar.uvmax, &stats))
    return -1;
/*
 * Get unit conversion factors, where needed.
 */
  switch(type->id) {
  case VS_PHS: /* Convert phases from radians to degrees */
    cnvfac = rtod;
    break;
  case VS_UMAG: /* Convert UV coordinates from wavelengths to the */
  case VS_VMAG: /*  currently selected units. */
  case VS_UVRAD:
    cnvfac = wavtouv(1.0);
    break;
  default:
    cnvfac = 1.0;
    break;
  };
/*
 * If this procedure is being called as a function, return the
 * statistics as an array of values.
 */
  if(outvals) {
/*
 * Allocate memory for the return array.
 */
    if((VOIDPTR(outvals)=valof_alloc(6,'f')) == NULL)
      return -1;
    outvals->adim[0] = 6;
    outvals->adim[1] = 1;
    outvals->adim[2] = 1;
    outvals->num_el = 6;
/*
 * Get a pointer to the return array, and copy the results into the array.
 */
    fptr = FLTPTR(outvals);
    fptr[0] = stats.nvis;
    fptr[1] = stats.mean * cnvfac;
    fptr[2] = stats.sigma * cnvfac;
    fptr[3] = stats.scatter * cnvfac;
    fptr[4] = stats.minval * cnvfac;
    fptr[5] = stats.maxval * cnvfac;
/*
 * If this procedure is being called as a command, print the results.
 */
  } else {
    lprintf(stdout, " N=%d Mean=%g +/- %g Scatter=%g Min=%g Max=%g\n",
	    stats.nvis, stats.mean * cnvfac, stats.sigma * cnvfac,
	    stats.scatter * cnvfac, stats.minval * cnvfac,
	    stats.maxval * cnvfac);
  };
  return no_error;
}

/*.......................................................................
 * Return the Rayleigh Jeans temperature that corresponds to a given flux
 * reading for a given planet at a given frequency and time.
 *
 * Input:
 *  flux      float     The flux measurement for the planet, at the
 *                      current mean frequency.
 *  planet     char *   Optional planet name to be used instead of the
 *                      name of the source of the current observation.
 *  epoch     float     The MJD UTC for which the flux was measured. If
 *                      omitted, the mean date of the observation is used.
 *  freq      float     The frequency at which the flux was measured. If
 *                      omitted, the mean frequency of the observation is
 *                      substituted.
 * Output:
 *  return    float     Optional return value gives the requested temperature.
 */
static Template(planet_temp_fn)
{
  float flux;           /* The flux to be converted */
  char *planet = NULL;  /* The name of the planet to be looked up */
  double mjd = 0.0;     /* The MJD UTC at which the flux was measured */
  double freq = 0.0;    /* The frequency to which the flux corresponds */
  double ra,dec;        /* The current geocentric RA,Dec of the planet */
  float diameter;       /* The angular diameter of the planet (radians) */
  float flattening;     /* The flattening of the planet (a-b)/a */
  float omega;          /* The solid angle subtended by the planet */
  float lambda;         /* The wavelength corresponding to 'freq' */
  float temp;           /* The brightness temperature */
/*
 * Get the parameters. Note that fall-throughs between cases.
 */
  switch(npar) {
  case 4:
    freq = *FLTPTR(invals[3]);
  case 3:
    mjd = *FLTPTR(invals[2]);
  case 2:
    planet = *STRPTR(invals[1]);
  case 1:
    flux = *FLTPTR(invals[0]);
    break;
  default:
    lprintf(stderr, "Wrong number of arguments.\n");
    return -1;
  };
/*
 * Fill in missing data from the current observation.
 */
  if((!planet || planet[0]=='\0' || freq <= 0.0 || mjd <= 0.0)) {
/*
 * Make sure that we have an observation to take data from.
 */
    if(nodata("planet_temp", OB_SELECT))
      return -1;
/*
 * Fill in the frequency?
 */
    if(freq <= 0.0)
      freq = getfreq(vlbob, -1);
/*
 * If no date was provided, substitute the mid time of the observation.
 */
    if(mjd <= 0.0) {
      double mid_ut = vlbob->rec[0].integ->ut +
	(vlbob->rec[vlbob->nrec-1].integ->ut - vlbob->rec[0].integ->ut)/2;
      long jd;
      double jdfrc, je;
      julday(mid_ut, vlbob->date.year, &jd, &jdfrc, &je);
      mjd = (jd - 2400000.5) + jdfrc;
    };
/*
 * If no source name was provided, use the source that is currently being
 * observed.
 */
    if(!planet || planet[0] == '\0')
      planet = vlbob->source.name;
  };
/*
 * Lookup the details of the source at the specified time and frequency.
 */
  if(planet_geometry(planet, mjd, &ra, &dec, &diameter, &flattening))
    return -1;
/*
 * Work out the solid angle subtended by the planet.
 */
  omega = pi * diameter * diameter * (1.0-flattening) / 4.0;
/*
 * Work out the wavelength that corresponds to freq.
 */
  lambda = cvel / freq;
/*
 * Compute the brightness temperature of the planet.
 */
  temp = flux/omega * 1.0e-26 * lambda * lambda / 2.0 / boltzmann;
/*
 * When called as a function, return the temperature. When called as
 * a command, report the temperature.
 */
  if(outvals)
    *FLTPTR(outvals) = temp;
  else
    lprintf(stdout, "Apparent brightness temperature = %g Kelvin\n", temp);
  return no_error;
}

/*.......................................................................
 * Return the geometry of a given planet at a given time.
 *
 * Input:
 *  planet     char *   Optional planet name to be used instead of the
 *                      name of the source of the current observation.
 *  epoch     float     The MJD UTC for which the flux was measured. If
 *                      omitted, the mean date of the observation is used.
 * Output:
 *  return    float     Optional return array gives the requested values.
 */
static Template(planet_geom_fn)
{
  char *planet = NULL;  /* The name of the planet to be looked up */
  double mjd = 0.0;     /* The MJD UTC at which the flux was measured */
  double ra,dec;        /* The current geocentric RA,Dec of the planet */
  float diameter;       /* The angular diameter of the planet (radians) */
  float omega;          /* The solid angle subtended by the planet */
  float flattening;     /* The flattening of the planet (a-b)/a */
/*
 * Get the parameters. Note that fall-throughs between cases.
 */
  switch(npar) {
  case 2:
    mjd = *FLTPTR(invals[1]);
  case 1:
    planet = *STRPTR(invals[0]);
    break;
  default:
    lprintf(stderr, "Wrong number of arguments.\n");
    return -1;
  };
/*
 * Fill in missing data from the current observation.
 */
  if((!planet || planet[0]=='\0' || mjd <= 0.0)) {
/*
 * Make sure that we have an observation to take data from.
 */
    if(nodata("planet_geometry", OB_SELECT))
      return -1;
/*
 * If no date was provided, substitute the mid time of the observation.
 */
    if(mjd <= 0.0) {
      double mid_ut = vlbob->rec[0].integ->ut +
	(vlbob->rec[vlbob->nrec-1].integ->ut - vlbob->rec[0].integ->ut)/2;
      long jd;
      double jdfrc, je;
      julday(mid_ut, vlbob->date.year, &jd, &jdfrc, &je);
      mjd = (jd - 2400000.5) + jdfrc;
    };
/*
 * If no source name was provided, use the source that is currently being
 * observed.
 */
    if(!planet || planet[0] == '\0')
      planet = vlbob->source.name;
  };
/*
 * Lookup the details of the source at the specified time and frequency.
 */
  if(planet_geometry(planet, mjd, &ra, &dec, &diameter, &flattening))
    return -1;
/*
 * Work out the solid angle subtended by the planet.
 */
  omega = pi * diameter * diameter * (1.0-flattening) / 4.0;
/*
 * If this procedure is being called as a function, return the
 * values as an array.
 */
  if(outvals) {
    float *fptr;     /* The return array */
/*
 * Allocate memory for the return array.
 */
    const int dim = 4;
    if((VOIDPTR(outvals)=valof_alloc(dim,'f')) == NULL)
      return -1;
    outvals->adim[0] = dim;
    outvals->adim[1] = 1;
    outvals->adim[2] = 1;
    outvals->num_el = dim;
/*
 * Get a pointer to the return array, and copy the results into the array.
 */
    fptr = FLTPTR(outvals);
    fptr[0] = diameter * rtoas;
    fptr[1] = flattening;
    fptr[2] = diameter * sqrt(1.0 - flattening) * rtoas;
    fptr[3] = omega * rtoas * rtoas;
/*
 * If this procedure is being called as a command, print the results.
 */
  } else {
    lprintf(stdout, "%s: Equatorial diameter=%g arcsec\n",
	    planet, diameter * rtoas, flattening);
    lprintf(stdout, "%*s  Flattening=%g\n", (int)strlen(planet), "",
	    flattening);
    lprintf(stdout, "%*s  Geometric diameter=%g arcsec\n",
	    (int)strlen(planet), "", diameter * sqrt(1.0 - flattening) * rtoas);
    lprintf(stdout, "%*s  Solid angle=%g arcsec^2\n", (int)strlen(planet), "",
	    omega * rtoas * rtoas);
  };
  return no_error;
}

/*.......................................................................
 * Return the Modified Julian Date that corresponds to a given Gregorian
 * date and time.
 *
 * Input:
 *  string    char *   A string containing a date and time, formatted
 *                     like dd-mmm-yyyy:hh:mm:ss.ss. Trailing time
 *                     components can be omitted.
 */
static Template(mjd_fn)
{
  double mjd;
/*
 * Attempt to parse the string.
 */
  if(parse_mjd(*STRPTR(invals[0]), 1, NULL, &mjd))
    return -1;
/*
 * Record the Modified Julian Date for return.
 */
  *FLTPTR(outvals) = mjd;
  return no_error;
}

/*.......................................................................
 * Replace the antenna voltage beams of a selection of antennas.
 *
 * Input:
 *  antennas    char * A list of antenna specifications.
 *  samples[]  float   The array of samples that define the antenna beam.
 *  binwidth   float   The radial width of each sample in samples[].
 *  freq       float   The frequency at which the primary beam is defined.
 */
static Template(antenna_beam_fn)
{
  int nsample = 0;      /* The number of samples */
  float binwidth = 0.0; /* The radial width covered by each sample */
  float *samples=NULL;  /* The array of samples */
  float freq = 0.0;     /* The frequency of the primary beam */
  char *spec=NULL;      /* The list of antenna specifications */
/*
 * Check that data have been read in.
 */
  if(nodata("antenna_beam", OB_INDEX))
    return -1;
/*
 * Get the parameters.
 */
  switch(npar) {
  case 4:
    freq = *FLTPTR(invals[3]);
    binwidth = xytorad(*FLTPTR(invals[2]));
    samples = FLTPTR(invals[1]);
    nsample = invals[1]->adim[0];
  case 1:
    spec = *STRPTR(invals[0]);
    break;
  default:
    lprintf(stderr, "Unexpected number of arguments.\n");
    return -1;
  };
/*
 * Assign the specified beam to the specified list of antennas.
 */
  if(set_antenna_beam(vlbob, spec, samples, nsample, binwidth, freq))
    return -1;
/*
 * Report what has been done.
 */
  if(npar==0)
    lprintf(stdout, "All antenna primary beams have been removed.\n");
  else
    lprintf(stdout, "A new antenna beam has been successfully installed.\n");
  return no_error;
}

/*.......................................................................
 * Set the pointing center of the observation.
 *
 * Input:
 *  ra,dec     The Right Ascension and Declination (sexagesimal strings).
 */
static Template(pointing_center_fn)
{
  char *ra_s = NULL;   /* The RA as a sexagesimal string */
  char *dec_s = NULL;  /* The Dec as a sexagesimal string */
  double ra, dec;      /* The RA and Dec in radians */
/*
 * Check that there is an observation to associate the position with.
 */
  if(nodata("pointing_center", OB_INDEX))
    return -1;
/*
 * Get the new Right Ascension and Declination.
 */
  switch(npar) {
  case 2:
    dec_s = *STRPTR(invals[1]);
    ra_s = *STRPTR(invals[0]);
/*
 * Parse the Right Ascension argument.
 */
    if(parse_sexagesimal_string(ra_s, &ra, NULL)) {
      lprintf(stderr, "pointing_center: Bad Right Ascension: %s\n", ra_s);
      return -1;
    };
    ra *= htor;
/*
 * Parse the Declination argument.
 */
    if(parse_sexagesimal_string(dec_s, &dec, NULL)) {
      lprintf(stderr, "pointing_center: Bad Declination: %s\n", dec_s);
      return -1;
    };
    dec *= dtor;
/*
 * Record the new values.
 */
    if(set_obs_radec(vlbob, ra, dec))
      return -1;
    break;
  case 0:   /* Lookup the existing values */
    ra = vlbob->source.obsra;
    dec = vlbob->source.obsdec;
    break;
  default:
    lprintf(stderr, "pointing_center: Missing declination.\n");
    return -1;
  };
/*
 * Report the current settings.
 */
  if(!vlbob->source.have_obs) {
    lprintf(stdout, "No pointing center is currently specified.\n");
  } else {
    char rabuf[15], decbuf[15];
    lprintf(stdout, "Pointing center:  RA=%s  Dec=%s  (%.1f)\n",
	    sradhms(ra, 3, 1, rabuf), sraddms(dec, 3, 1, decbuf),
	    vlbob->source.epoch);
  };
  return no_error;
}

/*.......................................................................
 * Replace the beams of all antennas to the square root of the specified
 * primary beam.
 *
 * Input:
 *  samples[]  float   The array of samples that define the primary beam.
 *  binwidth   float   The radial width of each sample in samples[].
 *  freq       float   The frequency at which the primary beam is defined.
 */
static Template(primary_beam_fn)
{
  int nsample = 0;      /* The number of samples */
  float binwidth = 0.0; /* The radial width covered by each sample */
  float *samples=NULL;  /* The array of samples */
  float freq = 0.0;     /* The frequency of the primary beam */
/*
 * Check that data have been read in.
 */
  if(nodata("primary_beam", OB_INDEX))
    return -1;
/*
 * Get the parameters.
 */
  switch(npar) {
  case 3:
    freq = *FLTPTR(invals[2]);
    binwidth = xytorad(*FLTPTR(invals[1]));
    samples = FLTPTR(invals[0]);
    nsample = invals[0]->adim[0];
  case 0:                         /* Remove the current primary beam */
    break;
  default:
    lprintf(stderr, "Unexpected number of arguments.\n");
    return -1;
  };
/*
 * Assign the specified beam to the specified list of antennas.
 */
  if(set_primary_beam(vlbob, samples, nsample, binwidth, freq))
    return -1;
/*
 * Report what has been done.
 */
  if(npar==0)
    lprintf(stdout, "All antenna primary beams have been removed.\n");
  else
    lprintf(stdout, "A new primary beam has been successfully installed.\n");
  return no_error;
}

/*.......................................................................
 * Flag the visibilites of one or more baselines.
 *
 * Input:
 *  spec        char *  A specification of a set of one or more baselines.
 *  doall    logical    If true, edit all channels/IFs. Otherwise only edit
 *                      the currently selected channels.
 *  ta,tb      float    The range of times to edit. If omitted, the
 *                      ta is given the start time of the data, and
 *                      tb the end time.
 */
static Template(flag_fn)
{
  char *spec="";       /* The telescope specification */
  int doall=0;         /* True to edit all channels/IFs */
  double mjd1=0,mjd2=0;/* The range of times to edit (MJD UTC) */
/*
 * Check that data have been selected.
 */
  if(nodata("flag_fn", OB_SELECT))
    return -1;
/*
 * Get the arguments.
 */
  switch(npar) {
  case 4:
    if(*STRPTR(invals[3])[0] != '\0' &&
       parse_mjd(*STRPTR(invals[3]), 1, NULL, &mjd2))
      return -1;
  case 3:
    if(*STRPTR(invals[2])[0] != '\0' &&
       parse_mjd(*STRPTR(invals[2]), 1, NULL, &mjd1))
      return -1;
  case 2:
    doall = *LOGPTR(invals[1]);
  case 1:
    spec = *STRPTR(invals[0]);
  };
/*
 * Perform the edits.
 */
  if(edit_baselines(vlbob, 1, spec, doall, mjd1, mjd2))
    return -1;
  return no_error;
}

/*.......................................................................
 * Un-flag the flagged visibilites of one or more baselines.
 *
 * Input:
 *  spec        char *  A specification of a set of one or more baselines.
 *  doall    logical    If true, edit all channels/IFs. Otherwise only edit
 *                      the currently selected channels.
 *  ta,tb      float    The range of times to edit. If omitted, the
 *                      ta is given the start time of the data, and
 *                      tb the end time.
 */
static Template(unflag_fn)
{
  char *spec="";        /* The telescope specification */
  int doall=0;          /* True to edit all channels/IFs */
  double mjd1=0,mjd2=0; /* The range of times to edit (MJD UTC) */
/*
 * Check that data have been selected.
 */
  if(nodata("unflag_fn", OB_SELECT))
    return -1;
/*
 * Get the arguments.
 */
  switch(npar) {
  case 4:
    if(*STRPTR(invals[3])[0] != '\0' &&
       parse_mjd(*STRPTR(invals[3]), 1, NULL, &mjd2))
      return -1;
  case 3:
    if(*STRPTR(invals[2])[0] != '\0' &&
       parse_mjd(*STRPTR(invals[2]), 1, NULL, &mjd1))
      return -1;
  case 2:
    doall = *LOGPTR(invals[1]);
  case 1:
    spec = *STRPTR(invals[0]);
  };
/*
 * Perform the edits.
 */
  if(edit_baselines(vlbob, 0, spec, doall, mjd1, mjd2))
    return -1;
  return no_error;
}

/*.......................................................................
 * Convert a distance in the map plane from the units that are currently
 * selected for the entry and display of map distances, to radians.
 *
 * Input:
 *  xy      float   The distance in map units.
 * Output:
 *  return  float   The distance in radians.
 */
static Template(map_to_rad_fn)
{
  *FLTPTR(outvals) = xytorad(*FLTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Convert a distance in the map plane from radians to the units that
 * are currently selected for the entry and display of map distances.
 *
 * Input:
 *  xy      float   The distance in map units.
 * Output:
 *  return  float   The distance in radians.
 */
static Template(rad_to_map_fn)
{
  *FLTPTR(outvals) = radtoxy(*FLTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Convert a distance in the uv plane from the units that are currently
 * selected for the entry and display of uv distances, to wavelengths.
 *
 * Input:
 *  uv      float   The distance in uv units.
 * Output:
 *  return  float   The distance in wavelengths.
 */
static Template(uv_to_wav_fn)
{
  *FLTPTR(outvals) = uvtowav(*FLTPTR(invals[0]));
  return no_error;
}

/*.......................................................................
 * Convert a distance in the uv plane from wavelengths to the units that
 * are currently selected for the entry and display of uv distances.
 *
 * Input:
 *  uv      float   The distance in uv units.
 * Output:
 *  return  float   The distance in wavelengths.
 */
static Template(wav_to_uv_fn)
{
  *FLTPTR(outvals) = wavtouv(*FLTPTR(invals[0]));
  return no_error;
}
