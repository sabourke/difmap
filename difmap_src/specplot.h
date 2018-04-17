#ifndef specplot_h
#define specplot_h

/* Enumerate possible X-axis units */

typedef enum {
  SP_CHAN,      /* Channel indexes */
  SP_FREQ       /* Frequency */
} SpXunit;

/* Enumerate X-axis recognized smoothing functions */

typedef enum {   /* Let v=channel, hwhm=HWHM of smoothing fn, and x=v/hwhm */
  SM_NONE,       /* f(x)=1.0 (No smoothing) */
  SM_HANNING,    /* f(x)=2*sin(pi*x)/(2*pi*x*(1-x)*(1+x)) */
  SM_GAUSSIAN,   /* f(x)=exp(-ln(2)*x^2) */
  SM_BOXCAR,     /* f(x)=1.0 for x=-1..1, and 0 elsewhere */
  SM_SINC        /* f(x)=sin(c*x)/(c*x), with c=1.8954942670340 */
} SmType;

typedef struct {
  SpXunit xunit;   /* The units of the smoothing width */
  SmType type;     /* The type of smoothing function */
  float fwhm;      /* The full-width at half the smoothing function maximum */
} SpSmooth;

/* Enumerate potentially variable spectrum indexing keys. */

typedef enum {
  SP_BASE,      /* Baseline selection */
  SP_POL,       /* Polarization selection */
  SP_TIME,      /* Time selection */
  SP_UVR,       /* UV radius selection */
  SP_NKEY       /* This should always be the last enumerator */
} SpKey;

/* Enumerate averaging modes */

typedef enum {
  SP_VECTOR,    /* Vector average visibilities versus time */
  SP_SCALAR     /* Scalar average visibilities versus time */
} SpAvMode;

/* Enumerate baseline selection modes */

typedef enum {
  SP_SPLIT,     /* Extract individual baselines from the first baselist */
  SP_GROUP      /* Display each baselist as a group */
} SpBMode;

/*
 * Group UV radius selection parameters. All radii are measured as
 * in wavelengths.
 * This is first initialized by new_Specattr() via sp_set_uvr().
 */

typedef struct {
  float uvrlim;   /* The max available UV radius */
  float uvmin;    /* The minimum UV radius to show spectra for */
  float uvmax;    /* The maximum UV radius to show spectra for */
  float uvstep;   /* The step size to break the range uvmin -> uvmax into */
} SpUV;

/*
 * There are so many specplot driving parameters that we need to be able
 * to set them up before calling specplot(). The following structure
 * will be initialized to defaults by new_Specattr() and can then be
 * modified via sp_*() method functions.
 */
typedef struct {
  double stime;     /* Start time */
  double etime;     /* End time */
  double scan;      /* Scan duration if >=0, or scan interval if < 0 */
  float amin,amax;  /* The amplitude range to plot */
  float pmin,pmax;  /* The phase range to plot */
  int ca, cb;       /* The channel range of the plot */
  int nplot;        /* The number of plots per page */
  int doamp;        /* True to plot amplitude spectrum */
  int dophs;        /* True to plot phase spectrum */
  int docross;      /* If true, use a cross-hair cursor by default */
  int dojoin;       /* If true, draw lines between channel values */
  int dohist;       /* Line join style: true -> histogram, false -> vector */
  int dobars;       /* True to plot error bars */
  Pollist *pl;      /* List of polarizations */
  Bgrplist *bgl;    /* List of baseline selection lists */
  SpUV uvr;         /* UV radius iterator parameters */
  SpKey key[SP_NKEY]; /* Array of spectrum indexing keys */
  int nkey;         /* The initial number of variable elements of key[] */
  SpXunit xunit;    /* X-axis units */
  SpSmooth smooth;  /* Smoothing parameters */
  SpBMode bmode;    /* Baseline selection mode */
  SpAvMode avmode;  /* Visibility averagin mode */
  Enumtab *xtsym;   /* Symbol-table of X-axis type names */
  Enumtab *keysym;  /* Symbol-table of selection key names */
  Enumtab *smsym;   /* Symbol-table of smoothing function names */
  Enumtab *avsym;   /* Symbol-table of averaging mode names */
  Enumtab *bmsym;   /* Symbol-table of baseline selection mode names */
} Specattr;

Specattr *new_Specattr(Observation *ob);
Specattr *del_Specattr(Specattr *sa);

int sp_set_pol(Observation *ob, Specattr *sa, Pollist *pol);
int sp_set_bgl(Observation *ob, Specattr *sa, SpBMode bmode, Bgrplist *bgl);
int sp_set_times(Observation *ob, Specattr *sa, double stime, double etime, 
		 double scan);
int sp_set_axes(Observation *ob, Specattr *sa, int ca, int cb,
		float amin, float amax, float pmin, float pmax);
int sp_set_options(Specattr *sa, int nplot, SpXunit xunit, SpAvMode avmode);
int sp_set_order(Specattr *sa, SpKey *keys, int nkey);
int sp_set_smooth(Specattr *sa, SpXunit xunit, SmType smtype, float fwhm);
int sp_set_uvrange(Observation *ob, Specattr *sa, float uvmin, float uvmax,
		   float uvstep);
char *sp_set_flags(Specattr *sa, char *flags);

int specplot(Observation *ob, Specattr *sa, int docurs, int npage);

#endif

