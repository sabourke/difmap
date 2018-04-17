#ifndef obs_h
#define obs_h

/*-------------------------------------------------------------------*/
/* Define the data structure to hold uv data read from a UV FITS file. */
/*-------------------------------------------------------------------*/

typedef struct Observation Observation;

#include "recio.h"
#include "dpage.h"
#include "ifpage.h"
#include "uvpage.h"
#include "model.h"
#include "chlist.h"
#include "pb.h"

#ifndef obedit_h
struct Obedit;
struct Edint;
#endif

/*
 * Define a source descriptor.
 */

typedef struct {
  char name[32];   /* Source name */
  double epoch;    /* The epoch of ra and dec */
  double ra;       /* Epoch mean right ascension of source */
  double dec;      /* Epoch mean declination of source */
  double app_ra;   /* Right ascension at date of observation */
  double app_dec;  /* Declination at date of observation */
  double tot_flux; /* Estimate of total flux in source (Jy) */
  int have_obs;    /* True if obsra and obsdec, are have been */
                   /* explicitly specified. If false, they are */
                   /* assigned the values of the ra and dec members. */
  double obsra;    /* The RA of the pointing-center (radians) */
  double obsdec;   /* The DEC of the pointing-center (radians) */
  float east;      /* The eastward offset of the source position (ra,dec) */
                   /*  from the pointing center (obsra,obsdec) in radians */
  float north;     /* The northward offset of the source position (ra,dec) */
                   /*  from the pointing center (obsra,obsdec) in radians */
} Source;

/*
 * Define a station descriptor, usable both for ground stations and
 * orbitting stations.
 */

typedef struct {
  char name[16];             /* Telescope name */
  int antno;                 /* FITS AN table antenna number */
  int antfix;                /* Flag to tell selfcal not to change gain */
  float antwt;               /* Extra weight to apply in selfcal */
  enum {GROUND, ORBIT} type; /* Station type */
  union {
    struct {     /* Ground based telecope parameters */
      double x;
      double y;
      double z;
    } gnd;
    struct {    /* Satellite based telescope parameters */
      double semi_major;  /* Semi-major axis of orbit (meters) */
      double eccentricity;/* Eccentricity of orbit */
      double inclination; /* Inclination of orbit to equator (degrees) */
      double ra_ascending;/* Right ascension of ascending node (degrees) */
      double arg_perigee; /* The argument of perigee (degrees) */
      double mean_anomoly;/* The mean anomaly at the reference time (degrees) */
    } orb;
  } geo;
  VoltageBeam *vb;        /* The voltage beam of the antenna, or NULL if not */
                          /* known. */
} Station;

/*
 * Struct to hold telescope phase and amplitude corrections.
 */
typedef struct {
  float amp_cor;   /* Amplitude correction */
  float phs_cor;   /* Phase correction */
  int bad;         /* If true, the correction has been flagged as bad */
} Telcor;

/*
 * Struct to hold baseline phase and amplitude corrections.
 */
typedef struct {
  float amp_cor;   /* Amplitude correction */
  float phs_cor;   /* Phase correction */
} Bascor;

/*
 * Struct to hold the sum of visibility weights on a given baseline.
 */
typedef struct {
  float wtsum;     /* The sum of visibility weights */
} Baswt;

/*
 * Define a struct used to contain memory to be distributed between
 * baseline descriptors. For example, the arrays of nif baseline
 * corrections in each descriptor is distributed from a single array
 * of nbase * nif elements held in the following struct.
 */
typedef struct {
  int nif;        /* The number of IFs currently allocated for */
  int nbase;      /* The number of baselines currently allocated for */
  Bascor *bcor;   /* nbase * nif baseline corrections */
  Baswt *bwt;     /* nbase * nif baseline weights */
} Basmem;

/*
 * Define the struct that will hold the details of one baseline.
 */
typedef struct {
  int tel_a, tel_b; /* Telescopes as 0-relative indexes into Station array */
  double boff;      /* Baseline hour angle (radians) */
  double bxy;       /* Baseline XY distance (meters) */
  double bz;        /* Baseline Z distance (meters) */
  Bascor *bcor;     /* nif IF-specific time-invariant gain corrections */
  Baswt *bwt;       /* nif IF-specific baseline weights. You must call */
                    /*  update_baseline_weights() before using this member. */
} Baseline;

/*
 * Define structs to be used to store otherwise unused items that
 * were read from binary table versions of AIPS AN antenna tables.
 * These are stored only to be re-producible when writing output
 * FITS files.
 */

typedef struct {
  double stabxyz[3];    /* Station X,Y,Z coords wrt arrayx,arrayy,arrayz */
  double *orbparm;      /* numorb orbital parameters for mntsta==2 */
  double staxof;        /* Axis offset (meters) */
  double polaa;         /* Feed A position angle (degrees) */
  double polab;         /* Feed B position angle (degrees) */
  double *polcala;      /* an->nopcal feed A polarization cal parameters */
  double *polcalb;      /* an->nopcal feed B polarization cal parameters */
  int mntsta;           /* Mount type, 0=altaz, 1=equatorial, 2=orbiting */
  int nosta;            /* Station number */
  char poltya;          /* Polarization type of feed A 'R','L','X','Y' */
  char poltyb;          /* Polarization type of feed B 'R','L','X','Y' */
  char anname[16];      /* Antenna name */
} Bintel;

typedef struct {
  double arrayx,arrayy,arrayz;  /* Array center */
  double gstia0;                /* GST at time=0 on ref date (degrees) */
  double degpdy;                /* Earth rotation rate (deg/day) */
  double freq;                  /* Sub-array reference freq */
  double polarx,polary;         /* Polar X,Y position at rdate (meters) */
  double ut1utc;                /* UT1-UTC (sec) */
  double datutc;                /* Data time - UTC (sec) */
  double *calpar;               /* Array of 2*nopcal*ob->nstat cal parameters */
  double *orbpar;               /* Array of numorb*ob->nstat cal parameters */
  Bintel *bt;                   /* ob->nstat telescope parameters */
  int nopcal;                   /* Number of polarization parameters */
  int numorb;                   /* Number of orbital parameters (0 or 6) */
  char arrnam[32];              /* Array name */
  char poltype[16];             /* Polarization parameterization type */
  char timsys[4];               /* Time system 'IAT' or 'UTC' */
  char rdate[9];                /* Sub-array reference date DD/MM/YY */
} Binan;

/*
 * Define bitmask union flagging types.
 */
typedef enum {
  FLAG_DEL=1,    /* Visibility deleted */
  FLAG_BAD=2,    /* Visibility flagged */
  FLAG_TA=4,     /* Telescope tel_a flagged (see Baseline) */
  FLAG_TB=8      /* Telescope tel_b flagged (see Baseline) */
} Flagtype;

/*
 * Data structure for one visibility.
 */
typedef struct {
  float amp;       /* Visibility amplitude */
  float modamp;    /* Model amplitude */
  float phs;       /* Visibility phase (radians) */
  float modphs;    /* Model phase (radians) */
  float wt;        /* Visibility weight (1/variance) */
  float u;         /* U coordinate (light-seconds) */
  float v;         /* V coordinate (light-seconds) */
  float w;         /* W coordinate (light-seconds) */
  float dt;        /* The integration time, or 0.0 if not known */
  int bad;         /* Status recorded as a union of Flagtype enumerators */
} Visibility;

/*
 * Define a container of IF-specific integration telescope corrections.
 */
typedef struct {
  Telcor *tcor;    /* Array of nstat telescope gain corrections */
} Intcor;

/*
 * Define a struct used to record blocks of memory that are distributed
 * between integration elements. For instance, rather than allocating
 * ntime separate arrays of nbase visibilities (one array per integration),
 * one array of nbase * ntime visibilities is allocated and distributed
 * between ntime integrations. This is more efficient in memory, helps
 * to prevent memory fragmentation and helps to prevent paging, by
 * ensuring that all visibilities come from the same part of virtual memory.
 */
typedef struct {
  int ntime;        /* The number of integrations catered for */
  int nbase;        /* The number of baselines catered for */
  int nstat;        /* The number of stations catered for */
  int nif;          /* The number of IFs catered for */
  Visibility *vis;  /* sub->ntime * sub->nbase visibilities */
  Intcor *icor;     /* sub->ntime * sub->nif Telcor array containers */
  Telcor *tcor;     /* sub->ntime * sub->nif * sub->nstat station corrections */
} Intmem;

/*
 * Define a struct to hold one integration.
 */

struct Subarray;

typedef struct Integration {
  double ut;            /* Start time of integration (seconds into year) */
  long irec;            /* Sequential record number in scratch files */
  struct Subarray *sub; /* Pointer to the associated sub-array descriptor */
  Visibility *vis;      /* Array of visibilities for this integration */
  Intcor *icor;         /* Array of nif arrays of nstat telescope corrections */
  struct Edint *edlist; /* List of edits to be applied to scratch files */
} Integration;

/*
 * Define a struct for encapsulating sub-array information.
 */

typedef struct Subarray {
  double scangap;     /* The gap used to delimit scans (seconds) */
  double datutc;      /* Original UV FITS time = UTC time + datutc */
  int nif;            /* Number of IFs == ob->nif */
  int ntime;          /* The number of integrations in this sub-array */
  int nstat;          /* The number of stations in this sub-array */
  int nbase;          /* The number of used baselines in this sub-array */
  Station *tel;       /* Array of 'nstat' station descriptors */
  Baseline *base;     /* Array of 'nbase' baseline descriptors */
  Binan *binan;       /* FITS binary AN-table sub-array info, or NULL */
  int p_refant;       /* 0-relative reference antenna number or -1 if not set */
  double *p_diff;     /* An array of 'nif' R-L phase differences for p_refant */
  Integration *integ; /* Array of sub->ntime integrations */
  struct Observation *ob;/* Pointer to parent Observation */
/*-----------------------------------------------------------------------*/
/* The following members are only for use by memory management functions */
/*-----------------------------------------------------------------------*/
  Intmem *imem;       /* Contains memory distributed between integrations */
  Basmem *bmem;       /* Contains memory distributed between baselines */
} Subarray;

typedef struct {
  double freq;    /* Frequency of first spectral-line channel (Hz) */
  double df;      /* Signed frequency offset between channels (Hz) */
  double bw;      /* Total bandwidth (Hz) */
  int coff;       /* Offset of channel 0 wrt set of all channels in all IFs */
  Chlist *cl;     /* List of channels within range 0..ob->nchan-1 */
                  /* If cl==NULL no channels are selected */
  int wtsum_bad;  /* This says that the visibility weights upon which */
                  /*  the per-baseline sums sub[*].base[*].wtsum[cif] */
                  /*  depend have changed in this IF since the sums were */
                  /*  last updated. This flag is set by calling */
                  /*  flag_baseline_weights() and used by */
                  /*  update_baseline_weights(). */
} If;

/*
 * Enumerate known stokes parameters and polarizations.
 */
typedef enum {
  NO_POL=0,
  SI=1, SQ=2, SU=3, SV=4,
  RR = -1, LL = -2, RL = -3, LR = -4, XX = -5, YY = -6, XY = -7, YX = -8,
  PI_POL = -9
} Stokes;

/* Stream polarization descriptor */
#define GETPOL_FN(fn) void (fn)(Obpol *pol, Cvis *pvis, Cvis *out)
typedef struct Obpol Obpol;
struct Obpol{
  Stokes type;   /* Current polarization/Stokes-parameter type. */
  int pa;        /* First of up to 2 indexes of stokes parameters to combine */
  int pb;        /* The second stokes parameter index, or -1 if not required */
  GETPOL_FN(*getpol); /* Use to extract a visibility of polarization 'type'*/
};

int get_Obpol(Observation *ob, Stokes stokes, int report, Obpol *obpol);

/*
 * Enumerate recognized spherical coordinate projections.
 */
typedef enum {
  PRJ_NON,      /* Unknown projection type */
  PRJ_SIN,      /* Sin (orthographic) projection */
  PRJ_NCP       /* North-Celestial-Pole (orthographic) projection. */
} Proj;

char *Proj_name(Proj proj);
Proj name_Proj(char *name);

/*
 * Declare a container to record parameters associated with the
 * specific selection of UV data used to construct the current
 * processing stream. 
 */
typedef struct {
  Chlist *cl;    /* List of combined spectral-line channel ranges */
  Obpol pol;     /* The descriptor of the stream polarization */
  int cif;       /* The index of the IF currently in memory */
  float uvscale; /* UVW coord multiplier for IF cif if ob->state>=OB_GETIF */
} UVstream;

/*
 * Declare a container to record an optional zero-spacing flux.
 */
typedef struct {
  float amp;      /* Estimated zero-spacing amplitude */
  float modamp;   /* Zero-spacing UV-model amplitude */
  float wt;       /* Visibility weight to assign the zero-spacing flux */
} UVzero;

/*
 * Declare a container to record any geometrical transformations
 * of the UV data in memory with respect to the data in the scratch files.
 */
typedef struct {
  float east;    /* Eastward shift applied to phases (radians) */
  float north;   /* Northward shift applied to phases (radians) */
  float uvangle; /* Clockwise rotation of UVW coordinates (radians) */
  float wtscale; /* Scale factor to be applied to weights */
} UVgeom;

/*
 * Define a container to hold the values of miscellaneous descriptive
 * keywords from a UV FITS header. String valued keywords that were not
 * provided in the input FITS file will be represented by a NULL pointer.
 */
typedef struct {
  char *origin;   /* Organization that created the FITS file */
  char *date_obs; /* Date of observation as DD/MM/YY */
  char *telescop; /* Telescope used in taking the observation */
  char *instrume; /* Instrument used in taking the observation */
  char *observer; /* The name of the observer */
  char *bunit;    /* Units of data */
  double equinox; /* The equinox of the data (0.0 if not in input FITS file) */
} Obhead;

/* Define a container for recording alternate velocity info */

typedef struct {   /* AIPS velocity information */
  int velref;      /* Velocity type 0 NONE >256 RADIO, 1 LSR, 2 HEL, 3 OBS */
  double altrval;  /* Reference velocity at pixel 'rpix' */
  double altrpix;  /* Reference pixel */
  double restfreq; /* Rest frequency (Hz) */
} Obvel;

/*
 * Define a struct to collect reference date related items.
 * NB. All integration UTs are measured as the number of seconds
 * since the start of year 'year'.
 */
typedef struct {
  int year;        /* The year of observation eg. 1992 */
  double utc_ref;  /* UTC date at start of year 'year' */  
  double ut;       /* UT of observation as seconds into year 'year' */
  double app_st;   /* Apparent Sidereal time at 'ut' */
  double cav_tim;  /* Coherent averaging time (s) */
  double iav_tim;  /* Incoherent averaging time (s) */
} Obdate;

/* Define a structure to associate integrations with scratch file records */

typedef struct {
  Integration *integ;
} Intrec;

/*
 * Enumerate the possible initialization states of the observation descriptor.
 */
typedef enum {
  OB_BAD,    /* Descriptor is corrupt */
  OB_ALLOC,  /* Descriptor allocated but no data read */
  OB_DATA,   /* Raw data copied to uvdata scratch file - set by get_uvdata() */
  OB_INDEX,  /* Data read and indexed - set by ini_Intrec() */
  OB_SELECT, /* Data selected with ob_select() */
  OB_RAWIF,  /* Uncorrected visibilities in memory [set by iniIF()] */
  OB_GETIF   /* Calibrated visibilities in memory [set by iniIF()] */
} Obstate;

/*
 * The main structure holding all details of the experiment including
 * source details, the array of stations and data arrays.
 */

struct Observation {
  Obstate state;   /* Current state of descriptor */
  int nhist;       /* Current number of history lines */
  int nsub;        /* The number of telescope sub-arrays */
  int nrec;        /* The total number of integrations in all sub-arrays */
  int nif;         /* Number of IFs */
  int npol;        /* Number of polarizations or stokes parameters */
  int nchan;       /* Number of spectral-line channels per IF */
  int nbmax;       /* The max number of baselines in any subarray */
  int nctotal;     /* The total number of channels in all IFs */
  int hasmod;      /* True if model visibilities exist */
  int have_inttim; /* True if the visibility integration times are usable */
  Obdate date;     /* Reference date info. */
  Obhead misc;     /* Miscellaneous values of descriptive FITS head info */
  Obvel vel;       /* AIPS alternate velocity definitions */
  Proj proj;       /* Spherical coordinate projection */
  UVstream stream; /* The UV data ranges used to construct the current stream */
  UVgeom geom;     /* Geometrical transformations applied to UV data */
  UVzero uvzero;   /* Zero-spacing flux */
  Source source;   /* Source characteristics */
  Stokes *pols;    /* Array of npol stokes parameter or polarization types */
  If *ifs;         /* Array of IF descriptors */
  Subarray *sub;   /* Array of 'nsub' sub-array descriptors */
  Intrec *rec;     /* Array of 'nrec' file record indexing descriptors */
  Dpage *dp;       /* uvdata paging descriptor */
  IFpage *ip;      /* IF paging descriptor */
  Recio *his;      /* History paging descriptor */
  UVpage *uvp;     /* UV model paging descriptor */
  Model *model;    /* The component form of the established UV model */
  Model *newmod;   /* The tentative, un-established part of the model */
  Model *cmodel;   /* Established continuum model */
  Model *cnewmod;  /* Un-established continuum model */
  struct ModelTable *mtab; /* A table of models corresponding to different */
                           /*  selections. */
  struct Obedit *obed; /* Pointer to container of deferred editing free-list */
  AntennaBeams *ab; /* A container of physical antenna and baseline beams */
};

/* Observation constructor */

Observation *new_Observation(const char *file_name, double binwid, int scatter,
			     int keepant, Chlist *cl, Stokes stokes);

/* Observation memory (re-)allocator */

Observation *Obs_alloc(Observation *ob, int ntotal, int nbmax, int nsub,
		       int nif, int npol, int nchan);

/* Observation destructor */

Observation *del_Observation(Observation *ob);

/* Check the validity of an Observation descriptor and report problems */

int ob_ready(Observation *ob, Obstate state, const char *name);

/* Select a new UV data combination */

int ob_select(Observation *ob, int keep, Chlist *cl, Stokes stokes);

/* Swap in a given IF into the associated Observation */

int getIF(Observation *ob, int cif);

/*
 * Use the following pair of calls to bracket code that potentially
 * changes the current stream IF.
 */
int get_cif_state(Observation *ob);
int set_cif_state(Observation *ob, int cif);

/* Swap in the UV model of a given IF */

int getmodel(Observation *ob, int cif);

/* Swap out the UV model of a given IF */

int putmodel(Observation *ob, int cif);

/* Defered IF index iterator */

int nextIF(Observation *ob, int cif, int skip_empty, int step);

/* Append a line of history */

int add_hist(Observation *ob, char *hisrec);

/* Effectively clear all history information */

int clr_hist(Observation *ob);

/* Display the current history */

int showhist(Observation *ob, int dopage);

/* Get the UVW coordinate scale factor for a given IF */

float getuvscale(Observation *ob, int cif);

/* Get the center freq and bandwidth of a given IF */

double getfreq(Observation *ob, int cif);
double getbw(Observation *ob, int cif);

/* Function to incrementally change the weight scale factor */

int wtscale(Observation *ob, float scale);

int visflags(Observation *ob, Visibility *vis, int nbase,
	     float uvmin, float uvmax, int *usable);

/* Shift the phase center of an observation */

int obshift(Observation *ob, float east, float north);
int obunshift(Observation *ob);

/* Shift the phase center of the observed visibilities only */

int uvshift(Observation *ob, float east, float north);

void uvrotate(Observation *ob, float angle);

typedef struct { /* Return value of uvrange() for the range uvmin -> uvmax */
  float uvrmin;  /* The min visibility UV radius (wavelengths) */
  float uvrmax;  /* The max visibility UV radius (wavelengths) */
  float umin;    /* The min visibility |U| distance (wavelengths) */
  float umax;    /* The max visibility |U| distance (wavelengths) */
  float vmin;    /* The min visibility |V| distance (wavelengths) */
  float vmax;    /* The max visibility |V| distance (wavelengths) */
  float ampmin;  /* The min visibility amplitude */
  float ampmax;  /* The max visibility amplitude */
  float wtmin;   /* The minimum visibility weight */
  float wtmax;   /* The maximum visibility weight */
} UVrange;

UVrange *uvrange(Observation *ob, int doall, int dores,
		 float uvmin, float uvmax);

/* Return a string of telescope abbreviations or array names */

int stnstr(Observation *ob, char *abrstr, int slen);

/* Clear the established model */

int clrmod(Observation *ob, int doold, int donew, int docont);

/* Add to the established model */

int obaddmod(Observation *ob, Model *mod, int keep, int docont, int append);

/* Establish the tentative model */

int mergemod(Observation *ob, int doold);

/* Append normal models to the continuum models */

int setcmod(Observation *ob, int tocont);

/* Partition the models into fixed established and variable tentative models */

int obvarmod(Observation *ob);

/* Edit the variable part or all of the established and tentative models */

int obedmod(Observation *ob, int dovar);

/* Add/subtract lone model components from the established or tentative model */

Modcmp *obaddcmp(Observation *ob, Modcmp *cmp, int keep);
Modcmp *obremcmp(Observation *ob, Modcmp *cmp, int keep);

/* Fit the variable part of the current models to the visibility residuals */

int fituvmodel(Observation *ob, int niter, float uvmin, float uvmax);

void vlbhead(Observation *ob);

void uncalib(Observation *ob, int doamp, int dophs, int doflag, int doreset);

int app_Telcor(Observation *ob, int cif);

int ed_Telcor(Observation *ob, Subarray *sub, int cif, int ut, int itel, int doflag);
int adj_Telcor(Observation *ob, Subarray *sub, int cif, int ut, int itel,
	       float amp_cor, float phs_cor);
int clr_Telcor(Observation *ob, Subarray *sub, int cif, int ut, int itel);


typedef struct { /* Container used to return values from moddif() */
  long ndata;    /* The number of measurements used */
  float uvmin;   /* The minimum of the selected range of UV radii */
  float uvmax;   /* The maximum of the selected range of UV radii */
  float chisq;   /* The (unreduced) value of chi-squared */
  float rms;     /* The RMS deviation between model and data (Jy) */
} Moddif;

int moddif(Observation *ob, Moddif *md, float uvmin, float uvmax);

int resoff(Observation *ob, int doall, int base, int isub);
int clroff(Observation *ob, int doall, int doamp, int dophs);

int ed_integ(Observation *ob, Subarray *sub, int ut, int cif, int doflag,
	     int selbase, int selstat, int selchan, int selif, int index);
int ed_flush(Observation *ob);

Observation *uvaver(Observation *ob, float avtime, int scatter);

char *Stokes_name(Stokes pol);
Stokes Stokes_id(char *name);

Observation *uvf_read(const char *name, double binwid, int scatter,
		      int keepant);
int uvf_write(Observation *ob, const char *name, int doshift);

/* If descriptor memory management functions */

If *new_If(Observation *ob, int nif);
If *del_If(Observation *ob);

/* Resoff baseline correction method functions */

int ini_bcor(Observation *ob, int doall, int doamp, int dophs);
int app_bcor(Observation *ob, int cif);

/* Binary FITS AN table descriptor memory management functions */

Binan *new_Binan(Subarray *sub, int nstat, int nopcal, int numorb);
Binan *del_Binan(Subarray *sub);
int fix_Binan(Subarray *sub, int *t_keep);

/* Sub-array memory management functions */

Subarray *new_Subarray(Observation *ob, int nsub);
int ini_Subarray(Subarray *sub, int nif, int nbase, int nstat,
		       int ntime);
Subarray *del_Subarray(Observation *ob);
int sub_bad(Subarray *sub, const char *fn);

/* Intrec memory management and initialization functions */

Intrec *new_Intrec(Observation *ob, int nrec);
int ini_Intrec(Observation *ob);
Intrec *del_Intrec(Observation *ob);

/* Obhead memory management functions */

int ini_Obhead(Observation *ob, char *origin, char *date_obs, char *telescop,
	       char *instrume, char *observer, char *bunit, double equinox);
void clr_Obhead(Observation *ob);

/* ob->dp service routines */

int dp_cal(Observation *ob);
int dp_shift(Observation *ob);
int dp_getpol(Observation *ob, Obpol *obpol);

/* Enumerate integration time-stamp search operators */

typedef enum {UT_LT, UT_LE, UT_NR, UT_GE, UT_GT} UTfind;

/* Return the ob->rec[] array index of a matching integration time-stamp */

int ob_find_ut(Observation *ob, double ut, UTfind op);

/* Return the sub->integ[] array index of a matching integration time-stamp */

int sub_find_ut(Subarray *sub, double ut, UTfind op);

/* Convert from direction cosine offsets to RA and DEC */

double lmtora(double ra, double dec, double l, double m, Proj proj);
double lmtodec(double ra, double dec, double l, double m, Proj proj);

/* Compute the direction cosine distances between two RA,DEC positions */

double radec_to_l(double ref_ra, double ref_dec, double ra, double dec,
		  Proj proj);
double radec_to_m(double ref_ra, double ref_dec, double ra, double dec,
		  Proj proj);

/*
 * Replace the current model with any model previously saved for
 * the currently selected channel-range and polarization.
 */
int ob_install_select_model(Observation *ob);

/*
 * Record the current model in the model table, indexed by the
 * currently selected channel-range and polarization.
 */
int ob_record_select_model(Observation *ob);

/*
 * Assign a given voltage beam to a given selection of antennas.
 */
int set_antenna_beam(Observation *ob, char *spec, float *samples, int nsample,
		     float binwidth, float freq);

/*
 * Set all antenna voltage beams to the square root of a given primary beam.
 */
int set_primary_beam(Observation *ob, float *samples, int nsample,
		     float binwidth, float freq);

/*
 * To allow primary beam calculations in observations that don't contain
 * OBSRA and OBSDEC keywords, specify the Right Ascension and Declination
 * of the original pointing center of the observation.
 */
int set_obs_radec(Observation *ob, double obsra, double obsdec);

/*
 * Get the offset of a given sky-projected position to the pointing center.
 */
float calc_pointing_offset(Observation *ob, float x, float y);

/*
 * Calculate the primary beam scale factor at a given map-projection
 * offset relative to the current map center, for a given baseline and
 * frequency.
 */
float pb_bl_factor(Subarray *sub, int base, double freq, float radius);

/*
 * Compute the primary beam scale factor averaged over all subarrays,
 * baselines and IFs for a given map-projection offset wrt the current
 * map center. This is the amount by which the flux of a delta function
 * at map position (x,y) is scaled by the combined effects of the primary
 * beams of all baselines.
 */
int pb_scale_factor(Observation *ob, float radius, float *factor);

/*
 * Correct the flux of a delta model component for the effects of the
 * combined primary beams of all baselines in all IFs.
 */
int pb_correct_delta_cmp(Observation *ob, Modcmp *cmp);

/*
 * Mark the per-baseline sums of visibility weights in a given IF as
 * stale wrt to the individual visibility weights of that IF. To flag
 * the sums in all IFs, specify cif as -1.
 */
int flag_baseline_weights(Observation *ob, int cif);

/*
 * Recompute the per-baseline sums of visibility weights in either a
 * given IF cif, or by specifying cif as -1, in all IFs.
 */
int update_baseline_weights(Observation *ob, int cif);

/*
 * Edit a set of baselines over a given time range.
 */
int edit_baselines(Observation *ob, int doflag, char *spec, int doall,
		   double mjd1, double mjd2);

#endif
