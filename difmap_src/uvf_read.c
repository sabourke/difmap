#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "logio.h"
#include "libfits.h"
#include "obs.h"
#include "vlbutil.h"
#include "vlbconst.h"
#include "intlist.h"
#include "visaver.h"
#include "slalib.h"

/*
 * Declare an internal struct used to record decoded random-parameter
 * values.
 */
typedef struct {
  double uu,vv,ww;   /* U,V and W coordinates in meters */
  double date;       /* The integration time-stamp (TAI Modified JD [days]) */
  double inttim;     /* The integration time, or 0.0 if not available */
  int isub;          /* FITS sub-array number (0-relative) */
  int ta,tb;         /* FITS station numbers */
  int fqid;          /* Frequency ID */
} Parval;

/*
 * Declare and define a container for random group parameter indexes.
 */
static struct {
  int ready;     /* Set to 1 only when the indexes have been initialized */
  int uu1,uu2;   /* Indexes of up to two UU random parameters */
  int vv1,vv2;   /* Indexes of up to two VV random parameters */
  int ww1,ww2;   /* Indexes of up to two WW random parameters */
  int bas1,bas2; /* Indexes of up to two BASELINE random parameters */
  int dat1,dat2; /* Indexes of up to two DATE random parameters */
  int fq1,fq2;   /* Indexes of up to two FREQSEL random parameters */
  int dt1,dt2;   /* Indexes of up to two INTTIM random parameters */
} gp = {0};

/*
 * Declare and define a container for PHDU axis indexes.
 */
static struct {
  int ready;     /* Set to 1 only when the indexes have been initialized */
  int cpos,cinc; /* Index and increment of COMPLEX axis */
  int spos,sinc; /* Index and increment of STOKES axis */
  int fpos,finc; /* Index and increment of FREQ axis */
  int ipos,iinc; /* Index and increment of IF axis */
  int rpos,rinc; /* Index and increment of RA axis */
  int dpos,dinc; /* Index and increment of DEC axis */
} ax = {0};

/*
 * Declare a type used to mark used stations or baselines
 * and for map from their indexes in the original AN table to their
 * indexes in the output sub-array.
 */
typedef struct {
  short used;     /* Set to 1 if used, 0 if not used */
  short slot;     /* Index of output slot in sub-array */
} Anmap;

/*
 * Declare a type to associate AN tables with sub-array and table descriptors
 * and record other sub-array specific info.
 */

typedef struct {
  Subarray *sub;      /* Pointer to sub-array descriptor */
  Thdu *thdu;         /* Pointer to AN table descriptor */
  Integration *integ; /* Pointer to next un-initiliazed integration in sub */
  int *antrow;        /* Row number (0..nrow-1) of station numbers (0..nsmax) */
  Anmap *bmap;        /* Map usage of the nbmax input baselines */
  Anmap *smap;        /* Map usage of stations by row number (0..nrow-1) */
  int ntime;          /* Number of integrations in this sub-array */
  int nrow;           /* The number of rows in the antenna table */
  int nsmax;          /* The highest station number cited in the table */
  int nstat;          /* The number of stations that are used */
  int nbmax;          /* Max possible number of baselines: nstat*(nstat-1)/2 */
  int nbase;          /* Number of used baselines */
  double datutc;      /* The data time - UTC */
} Antab;

/* Define the structure used to store details about a random-groups
 * UV fits file.
 */
typedef struct Fitob {
  Fits *fits;         /* The descriptor of the UV FITS file */
  int npar;           /* Number of random-parameters per group */
  double *pars;       /* Buffer array of 'npar' elements */
  int ndata;          /* Number of elements per group array */
  double *data;       /* Buffer array of 'ndata' elements */
  int maxan;          /* Max AN table version number */
  Antab *antab;       /* [0..maxan] Map AN version to internal descriptors. */
  int nbmax;          /* Max number of baselines per sub-array */
  int nsub;           /* The number of telescope sub-arrays */
  int ntime;          /* The total number of integrations in all sub-arrays */
  int nif;            /* The number of IFs */
  int npol;           /* The number of polarizations */
  int nchan;          /* The number of spectral-line channels */
  int scatter;        /* If true substitute scatter estimates of weights */
  double binwid;      /* Integration bin width (seconds). */
  double wtsign;      /* The sign of the AIPS WTSCAL factor */
  Intlist *ilist;     /* Integration bin list/iterator */
  Obdate date;        /* Observation reference date info recorded by get_date */
  Proj proj;          /* The UU,VV,WW projection code */
  long start_group;   /* The index of the first group with a useable date */
} Fitob;

static Fitob *new_Fitob(const char *name, double binwid, int scatter,
			int keepant);
static Fitob *del_Fitob(Fitob *fob);

static int loc_par(char *name, Phdu *phdu, int need, int fixlen,
		   int *pa, int *pb);
static int loc_axis(char *name, Phdu *phdu, int need, int *ax, int *inc);

static int grp_parms(Fitob *fob);
static int uvw_parms(Fitob *fob);
static int get_axes(Fitob *fob);
static double *get_data(Fitob *fob, long group);
static int get_source(Observation *ob, Fits *fits, Obdate *date);
static int find_subarrays(Fitob *fob, int keepant);
static int count_IFs(Fitob *fob);
static int count_stokes(Fitob *fob);
static int count_sources(Fitob *fob);
static int count_FQ_entries(Fitob *fob);
static int get_date(Fitob *fob);
static int get_subarray_time_systems(Fitob *fob, double iatutc);
static int get_antrow(Fits *fits, Thdu *thdu, Antab *an);
static int get_stations(Observation *ob, Fits *fits, Antab *an);
static int getanbin(Observation *ob, Fits *fits, Antab *an);
static int getanasc(Observation *ob, Fits *fits, Antab *an);
static int rd_p_refant(Fits *fits, Hdu *hdu, Subarray *sub);
static Intlist *bin_uvdata(Fitob *fob, double binwid);

static int get_misc(Observation *ob, Fits *fits);
static int get_vel(Observation *ob, Fits *fits);
static int get_IF_freq(Observation *ob, Fits *fits);
static int get_stokes(Observation *ob, Fits *fits);
static int get_history(Observation *ob, Fits *fits, Fitob *fob);
static int check_history(Observation *ob, char *hline, Fitob *fob);
static int get_baselines(Observation *ob, Fitob *fob, Antab *an);
static int get_subarrays(Observation *ob, Fitob *fob);
static int get_uvdata(Observation *ob, Fitob *fob);
static Parval *read_pars(Fitob *fob, long group);
static Anmap *loc_base(Fitob *fob, Parval *pval);

static Observation *foberr(Fitob *fob, Observation *ob);
static int uvretfn(Visaver *av, int iret);
static int string_is_empty(char *string);

/*.......................................................................
 * Read a new observation from a FITS file.
 *
 * Input:
 *  name          char *  The name of a random-group UV FITS file.
 *  binwid      double    The integration bin width to collect visibilities
 *                        into (seconds). If binwid < 1.0 no binning will
 *                        be performed.
 *  scatter        int    If true substitute weights deduced from the
 *                        scatter of data within each integration bin,
 *                        for the read data weights.
 *  keepant        int    If true, allocate space for all antennas and
 *                        associated baselines. If false, discard all
 *                        antennas and baselines that don't have
 *                        visibilities associated with them.
 * Output:
 *  return Observation *  The pointer to the new descriptor, or NULL on
 *                        error.
 */
Observation *uvf_read(const char *name, double binwid, int scatter, int keepant)
{
  Observation *ob; /* The descriptor to be returned */
  Fitob *fob;      /* The FITS/Observation intermediary descriptor */
/*
 * Check arguments.
 */
  if(name==NULL) {
    lprintf(stderr, "uvf_read: NULL file name intercepted.\n");
    return NULL;
  };
/*
 * Open the FITS file and obtain sufficient information to allocate
 * an Observation structure.
 */
  fob = new_Fitob(name, binwid, scatter, keepant);
  if(fob==NULL)
    return NULL;
/*
 * Create an observation descriptor sufficient to contain a single
 * solution-time's worth of data.
 */
  ob = Obs_alloc(NULL, fob->ntime, fob->nbmax, fob->nsub, fob->nif, fob->npol,
		 fob->nchan);
  if(ob==NULL)
    return foberr(fob, ob);
/*
 * Record the miscellaneous descriptive header keyword values.
 */
  if(get_misc(ob, fob->fits))
    return foberr(fob, ob);
/*
 * Record AIPS altdef velocity info if given.
 */
  if(get_vel(ob, fob->fits))
    return foberr(fob, ob);
/*
 * Are integration times available?
 */
  ob->have_inttim = gp.dt1 >= 0;
/*
 * Record the reference date details in ob.
 */
  ob->date = fob->date;
/*
 * Record the spherical projection type of the UVW coordinates.
 */
  ob->proj = fob->proj;
/*
 * Determine and record source characteristics in ob->source.
 */
  if(get_source(ob, fob->fits, &fob->date))
    return foberr(fob, ob);
/*
 * Intialize sub-array descriptors from AN tables.
 */
  if(get_subarrays(ob, fob))
    return foberr(fob, ob);
/*
 * Get IF frequency info.
 */
  if(get_IF_freq(ob, fob->fits))
    return foberr(fob, ob);
/*
 * Get polarization info.
 */
  if(get_stokes(ob, fob->fits))
    return foberr(fob, ob);
/*
 * Store FITS history.
 */
  if(get_history(ob, fob->fits, fob))
    return foberr(fob, ob);
/*
 * Read the UV data.
 */
  if(get_uvdata(ob, fob))
    return foberr(fob, ob);
/*
 * The FITS file and its Fitob intermediary descriptor are no
 * longer required.
 */
  fob = del_Fitob(fob);
/*
 * Return the initialized Observation.
 */
  return ob;
}

/*.......................................................................
 * Private function of uvf_read(), used to clean up after errors.
 */
static Observation *foberr(Fitob *fob, Observation *ob)
{
  fob = del_Fitob(fob);
  return del_Observation(ob);
}

/*.......................................................................
 * Open and interpret the header and tables of a random-group UV FITS file.
 * Record the results in a Fitob descriptor to be subsequently used by
 * read_fits().
 *
 * Input:
 *  name          char *  The name of a random-group UV FITS file.
 *  binwid      double    The integration bin width to collect visibilities
 *                        into (seconds). If binwid < 1.0 no binning will
 *                        be performed.
 *  scatter        int    If true substitute weights deduced from the
 *                        scatter of data within each integration bin,
 *                        for the read data weights.
 *  keepant        int    If true, allocate space for all antennas and
 *                        associated baselines. If false, discard all
 *                        antennas and baselines that don't have
 *                        visibilities associated with them.
 * Output:
 *  return       Fitob *  The pointer to the new descriptor, or NULL on
 *                        error.
 */
static Fitob *new_Fitob(const char *name, double binwid, int scatter,
			int keepant)
{
  Fitob *fob;  /* The return descriptor */
  int *dims;   /* Pointer into FITS dims array */
  Fits *fits;  /* Pointer to the FITS descriptor */
  Phdu *phdu;  /* Pointer to the PRIMARY HDU descriptor */
  int nsource; /* The number of sources in the file */
  int nfq;     /* The number of frequency group definitions in the file */
  int i;
/*
 * Allocate and zero initialize the required Fitob descriptor.
 */
  fob = (Fitob *) malloc(sizeof(Fitob));
  if(fob==NULL) {
    lprintf(stderr, "new_Fitob: Insufficient memory for descriptor.\n");
    return NULL;
  };
/*
 * Clear the descriptor.
 */
  fob->ndata = 0;
  fob->npar = 0;
  fob->nbmax = 0;
  fob->nsub = 0;
  fob->ntime = 0;
  fob->nif = 0;
  fob->npol = 0;
  fob->nchan = 0;
  fob->maxan = 0;
  fob->data = NULL;
  fob->antab = NULL;
  fob->pars = NULL;
  fob->fits = NULL;
  fob->scatter = scatter;
  fob->binwid = binwid < 1.0 ? 0.0 : binwid;
  fob->wtsign = 1.0;
  fob->ilist = NULL;
  fob->proj = PRJ_SIN;
  fob->start_group = 0L;
/*
 * Attempt to open the new FITS file.
 */
  fits = fob->fits = new_Fits(name, 1, 1, 0, 1);
  if(fits==NULL)
    return del_Fitob(fob);
/*
 * Keep user informed.
 */
  lprintf(stdout, "Reading UV FITS file: %s\n", name);
/*
 * Get the descriptor of the PRIMARY HDU and its dimensions.
 */
  phdu = (Phdu *) fits->hdu;
  dims = phdu->dims;
/*
 * The primary HDU must be a random-group HDU.
 */
  if(!phdu->groups || phdu->pcount==0) {
    lprintf(stderr, "get_fits: Error: Primary header does not contain random-groups.\n");
    return del_Fitob(fob);
  };
/*
 * Get the 0-relative indexes of each of the recognized random parameters
 * and axes - also return the size of the data array in fob->ndata.
 */
  if(grp_parms(fob) || get_axes(fob))
    return del_Fitob(fob);
/*
 * Determine the size of the group data-array, and allocate a data array
 * of this size to read a single group-array into
 */
  fob->ndata = 1;
  for(i=1; i<phdu->naxis; i++)
    fob->ndata *= dims[i];
  fob->data = (double *) malloc(sizeof(double) * fob->ndata);
/*
 * Also allocate an array to read random-group parameters into.
 */
  fob->npar = phdu->pcount;
  fob->pars = (double *) malloc(sizeof(double) * fob->npar);
/*
 * Insufficient memory?
 */
  if(fob->pars==NULL || fob->data==NULL) {
    lprintf(stderr, "new_Fitob: Insufficient memory to read FITS file.\n");
    return del_Fitob(fob);
  };
/*
 * Count the number of sources in the file.
 */
  nsource = count_sources(fob);
  if(nsource <= 0)
    return del_Fitob(fob);
/*
 * Too many sources?
 */
  if(nsource>1) {
    lprintf(stderr, "Unable to handle multi-source files.\n");
    return del_Fitob(fob);
  };
/*
 * Count the number of frequency groups in the file.
 */
  nfq = count_FQ_entries(fob);
  if(nfq <= 0)
    return del_Fitob(fob);
/*
 * Too many frequency groups?
 */
  if(nfq>1) {
    lprintf(stderr, "Unable to handle multi-frequency files.\n");
    return del_Fitob(fob);
  };
/*
 * Determine the number of stations and sub-arrays involved in the
 * observation.
 */
  if(find_subarrays(fob, keepant))
    return del_Fitob(fob);
/*
 * Work out the number of IFs.
 */
  if(count_IFs(fob) <= 0)
    return del_Fitob(fob);
/*
 * Work out the number of polarizations or stokes parameters recorded.
 */
  if(count_stokes(fob) <= 0)
    return del_Fitob(fob);
/*
 * Determine the number of spectral-line channels per IF.
 */
  fob->nchan = dims[ax.fpos];
/*
 * Read the first group to determine the start date and use this to
 * fill in the reference date info in fob->date.
 */
  if(get_date(fob))
    return del_Fitob(fob);
/*
 * Read through the UV data to associate groups into integrations,
 * to count the number of such in each sub-array, to record the date
 * of the first integration, and for each AN table to record which tables
 * are used, which baselines and antennas are used, their output baseline
 * and station indexes and their numbers.
 */
  fob->ilist = bin_uvdata(fob, fob->binwid);
  if(fob->ilist==NULL)
    return del_Fitob(fob);
/*
 * Return the initialized descriptor.
 */
  return fob;
}

/*.......................................................................
 * Delete a Fitob descriptor.
 *
 * Input:
 *  fob     Fitob *  The descriptor to be deleted.
 * Output:
 *  return  Fitob *  Allways NULL.
 */
static Fitob *del_Fitob(Fitob *fob)
{
  if(fob==NULL)
    return fob;
  if(fob->pars)
    free(fob->pars);
/*
 * Delete temporary AN table indexing array.
 */
  if(fob->antab) {
    Antab *an = fob->antab;
    int i;
    for(i=0; i<fob->maxan; i++,an++) {
      if(an->antrow)
	free(an->antrow);
      if(an->bmap)
	free(an->bmap);
      if(an->smap)
	free(an->smap);
    };
    free(fob->antab);
  };
  if(fob->data)
    free(fob->data);
  if(fob->fits)
    fob->fits = del_Fits(fob->fits);
  free(fob);
  return NULL;
}

/*.......................................................................
 * Record miscellaneous parameters from the FITS primary header.
 *
 * Input:
 *  ob Observation *  The descriptor of the Observation being read.
 *  fits      Fits *  The descriptor of the FITS file.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int get_misc(Observation *ob, Fits *fits)
{
  Phdu *phdu = (Phdu *)fits->hdu;  /* Descriptor of primary HDU */
/*
 * Record the keyword values.
 */
  if(ini_Obhead(ob, phdu->origin, phdu->date_obs, phdu->telescop,
		phdu->instrume, phdu->observer, phdu->bunit, phdu->equinox))
    return 1;
  return 0;
}

/*.......................................................................
 * Record the AIPS ALTDEF velocity information from the header in ob->vel.
 *
 * Input:
 *  ob Observation *  The descriptor of the Observation being read.
 *  fits      Fits *  The descriptor of the FITS file.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int get_vel(Observation *ob, Fits *fits)
{
  Fitkey key;        /* Descriptor of FITS keyword/value pair */
/*
 * Create a table of known velocity keywords.
 */
  enum {VELREF, ALTRVAL, ALTRPIX, RESTFREQ};
  static Fitkey velkeys[]={
    {"VELREF",  0, VELREF,   DAT_INT, NULL, NULL},
    {"ALTRVAL", 0, ALTRVAL,  DAT_DBL, NULL, NULL},
    {"ALTRPIX", 0, ALTRPIX,  DAT_DBL, NULL, NULL},
    {"RESTFREQ",0, RESTFREQ, DAT_DBL, NULL, NULL},
  };
/*
 * Rewind the primary header.
 */
  new_hline(fits->hdu, 0);   /* Rewind header */
/*
 * Search for the optional velocity keywords and read their values
 * where found.
 */
  while(next_key(fits, fits->hdu, velkeys,
		 sizeof(velkeys)/sizeof(Fitkey), EOH_SEEK, &key) == 0) {
    switch(key.keyid) {
    case VELREF:
      ob->vel.velref = KEYINT(key);
      break;
    case ALTRVAL:
      ob->vel.altrval = KEYDBL(key);
      break;
    case ALTRPIX:
      ob->vel.altrpix = KEYDBL(key);
      break;
    case RESTFREQ:
      ob->vel.restfreq = KEYDBL(key);
      break;
    };
  };
  return 0;
}

/*.......................................................................
 * Get the positions of up to two versions of a random parameter.
 *
 * Input:
 *  name   char *    The name of the random group parameter to locate.
 *  phdu   Phdu *    The primary HDU descriptor.
 *  need    int      If true then an error message will be evoked if the
 *                   parameter is not found.
 *  fixlen  int      If > 0, then this defines the maximum number of characters
 *                   to be compared. This makes it possible to search by
 *                   prefixes.
 * Input/Output:
 *  pa      int *    If pa!=NULL *pa will be assigned the 0-relative index
 *                   of the first matching parameter, or -1 if not found.
 *  pb      int *    If pb!=NULL *pb will be assigned the 0-relative index
 *                   of the second matching parameter, or -1 if not found.
 * Output:
 *  return  int      0 - OK.
 *                   1 - Parameter not found.
 */
static int loc_par(char *name, Phdu *phdu, int need, int fixlen,
		   int *pa, int *pb)
{
  int ia=0;   /* 1-relative index of fist matching parameter */
  int ib=0;   /* 1-relative index of second matching parameter */
/*
 * Parameter not yet found.
 */
  if(pa) *pa = -1;
  if(pb) *pb = -1;
/*
 * Sanity check.
 */
  if(name==NULL || phdu==NULL) {
    lprintf(stderr, "loc_par: NULL %s intercepted\n",
	    name==NULL ? "name" : "Phdu");
    return 1;
  };
/*
 * Locate the first matching parameter.
 */
  ia = find_gpar(phdu, name, fixlen, 1);
/*
 * See if the parameter is cited twice.
 */
  if(ia>0)
    ib = find_gpar(phdu, gpar_name(phdu, ia), 0, ia+1);
/*
 * Match?
 */
  if(ia>0) {
    if(pa) *pa = ia-1;
    if(pb) *pb = ib-1;
    return 0;
  };
/*
 * No match.
 */
  if(need) {
    lprintf(stderr,
	    "loc_par: Unable to locate required %s random parameter.\n", name);
  };
  return 1;
}

/*.......................................................................
 * Get the positions of required axes in the primary HDU.
 *
 * Input:
 *  name   char *    The name of the random group parameter to locate.
 *  phdu   Phdu *    The primary HDU descriptor.
 *  need    int      If true then an error message will be evoked if the
 *                   parameter is not found.
 * Input/Output:
 *  ax      int *    If ax!=NULL *ax will be assigned the 0-relative index
 *                   of the first matching axis, or -1 if not found.
 *  inc     int *    If inc!=NULL *inc will be assinged The increment to
 *                   the next element on the axis. 
 * Output:
 *  return  int      0 - OK.
 *                   1 - Axis not found.
 */
static int loc_axis(char *name, Phdu *phdu, int need, int *ax, int *inc)
{
  int iax;   /* The 1-relative index of the axis */
/*
 * Axis not yet found.
 */
  if(ax) *ax = -1;
  if(inc) *inc = 0;
/*
 * Sanity check.
 */
  if(name==NULL || phdu==NULL) {
    lprintf(stderr, "loc_par: NULL %s intercepted\n",
	    name==NULL ? "name" : "Phdu");
    return 1;
  };
/*
 * Locate the first matching axis.
 */
  iax = find_axis(phdu, name, 0, 1);
/*
 * Match?
 */
  if(iax>0) {
    if(ax) *ax = iax-1;
/*
 * Work out the increment between elements on the given axis.
 */
    if(inc) {
      int iinc=1;  /* The increment to be returned */
      int *dims = phdu->dims;
      int i;
      for(i=1; i<iax-1; i++)
	iinc *= dims[i];
      *inc = iinc;
    };
    return 0;
  };
/*
 * No match.
 */
  if(need) {
    lprintf(stderr,
	    "loc_axis: Unable to locate required %s axis.\n", name);
  };
  return 1;
}

/*.......................................................................
 * Process the array axes of a random-groups UV FITS file.
 *
 *  fob    Fitob *  The FITS-Observation intermediary descriptor.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int get_axes(Fitob *fob)
{
  Phdu *phdu;   /* Pointer to the primary HDU descriptor. */
/*
 * Axis increments are not yet useable.
 */
  ax.ready = 0;
/*
 * Get the primary HDU descriptor.
 */
  phdu = (Phdu *) fob->fits->hdu;
/*
 * Locate the required axes.
 */
  if(loc_axis("COMPLEX", phdu, 1, &ax.cpos, &ax.cinc) ||
     loc_axis("FREQ", phdu, 1, &ax.fpos, &ax.finc) ||
     loc_axis("RA", phdu, 1, &ax.rpos, &ax.rinc) ||
     loc_axis("DEC", phdu, 1, &ax.dpos, &ax.dinc))
    return 1;
/*
 * Locate optional axes.
 */
  loc_axis("STOKES", phdu, 0, &ax.spos, &ax.sinc);
  loc_axis("IF", phdu, 0, &ax.ipos, &ax.iinc);
/*
 * The COMPLEX axis MUST be the first axis. (AIPS FITTP and UVLOD
 * require this).
 */
  if(ax.cpos != 1) {
    lprintf(stderr, "get_axes: Illegal CTYPE1 (should be COMPLEX).\n");
    return 1;
  };
/*
 * The COMPLEX axis must have 3 elements (real,imag,weight).
 */
  if(phdu->dims[1] != 3) {
    lprintf(stderr,
	    "get_axes: COMPLEX axis has %d elements, it ought to have 3.\n",
	    phdu->dims[1]);
    return 1;
  };
/*
 * Increments are now initialized for use.
 */
  ax.ready = 1;
  return 0;
}

/*.......................................................................
 * Read an SU table or the primary header to determine various source
 * characteristics. Record the results in ob->source.
 * This function also reads the source frequency offsets for each IF and
 * adds them to ob->if[IF].freq.
 *
 * Input:
 *  ob   Observation *  The descriptor being initialized.
 *  fits        Fits *  The descriptor of the FITS file.
 *  date      Obdate *  The reference date descriptor.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int get_source(Observation *ob, Fits *fits, Obdate *date)
{
  Source *sou;         /* Pointer into ob->source. */
  Phdu *phdu;          /* The primary HDU */
  Thdu *thdu;          /* The SU table descriptor. */
  int icol;            /* 1-relative column index */
  int nchar;           /* Number of characters read */
  int have_obsra=0;    /* True if the OBSRA keyword was found */
  int have_obsdec=0;   /* True if the OBSDEC keyword was found */
  double obsra=0.0;    /* The Right Ascension of the pointing center */
  double obsdec=0.0;   /* The Declination of the pointing center */
  Fitkey key;          /* The lastest keyword/value pair found */
  int i;
/*
 * List the optional primary HDU keywords that are to be looked for.
 */
  enum {OBSRA, OBSDEC};
  static Fitkey misc_keys[]={
    {"OBSRA",  0, OBSRA,   DAT_DBL, NULL, NULL},
    {"OBSDEC", 0, OBSDEC,  DAT_DBL, NULL, NULL},
  };
/*
 * List required SU-table fields.
 */
  enum {COL_SOURCE, COL_RAEPO, COL_DECEPO, COL_EPOCH, COL_RAAPP, COL_DECAPP,
	COL_IFLUX, COL_FREQOFF};
  typedef struct {
    char *name;    /* Name of table column */
    int icol;      /* 1-relative column number */
  } Col;
  Col cols[] = {
    {"SOURCE",  0},
    {"RAEPO",   0},
    {"DECEPO",  0},
    {"EPOCH",   0},
    {"RAAPP",   0},
    {"DECAPP",  0},
    {"IFLUX",   0},
    {"FREQOFF", 0},
  };
  Col *col;  /* Pointer into cols[] */
/*
 * Sanity checks.
 */
  if(ob==NULL || fits==NULL) {
    lprintf(stderr, "get_source: NULL %s descriptor intercepted.\n",
	    ob==NULL ? "Observation":"Fits");
    return 1;
  };
/*
 * Get a pointer to the PRIMARY HDU descriptor.
 */
  phdu = (Phdu *) fits->hdu;
/*
 * Get a pointer to the source descriptor.
 */
  sou = &ob->source;
/*
 * Look for an SU table.
 */
  thdu = find_table(fits, "AIPS SU", 0, NULL);
/*
 * If there is a table, make sure that it only contains a single source
 * then read the source parameters.
 */
  if(thdu!=NULL) {
/*
 * Trap multi-source files.
 */
    if(numrow(thdu) != 1) {
      lprintf(stderr, "get_source: Unable to handle multi-source files.\n");
      return 1;
    };
/*
 * Locate all the required table columns.
 */
    for(icol=0,col=cols; icol<sizeof(cols)/sizeof(Col); icol++,col++) {
      col->icol = find_column(thdu, col->name, 0);
      if(col->icol==0) {
	lprintf(stderr,
		"get_source: Failed to find %s column in AIPS SU table.\n",
		col->name);
	return 1;
      };
    };
/*
 * Get the source name.
 */
    nchar = rcolumn(fits, thdu, cols[COL_SOURCE].icol, 1, DAT_CHR, 1, NULL, 0,
		    (int) sizeof(sou->name), sou->name);
    if(nchar < 1)
      return 1;
    stripstr(sou->name, nchar);  /* Remove trailing spaces */
/*
 * Get the RA (epoch) of the source, converted to radians.
 */
    if(rcolumn(fits, thdu, cols[COL_RAEPO].icol, 1, DAT_DBL, 1, NULL, 0, 1,
	       &sou->ra) != 1)
      return 1;
    sou->ra *= dtor;
/*
 * Get the DEC (epoch) of the source, converted to radians.
 */
    if(rcolumn(fits, thdu, cols[COL_DECEPO].icol, 1, DAT_DBL, 1, NULL, 0, 1,
	       &sou->dec) != 1)
      return 1;
    sou->dec *= dtor;
/*
 * Get the EPOCH of RA and DEC of the source.
 */
    if(rcolumn(fits, thdu, cols[COL_EPOCH].icol, 1, DAT_DBL, 1, NULL, 0, 1,
	       &sou->epoch) != 1)
      return 1;
/*
 * Get the RA (apparent) of the source, converted to radians.
 */
    if(rcolumn(fits, thdu, cols[COL_RAAPP].icol, 1, DAT_DBL, 1, NULL, 0, 1,
	       &sou->app_ra) != 1)
      return 1;
    sou->app_ra *= dtor;
/*
 * Get the DEC (apparent) of the source, converted to radians.
 */
    if(rcolumn(fits, thdu, cols[COL_DECAPP].icol, 1, DAT_DBL, 1, NULL, 0, 1,
	       &sou->app_dec) != 1)
      return 1;
    sou->app_dec *= dtor;
/*
 * Get the total flux of the source.
 */
    if(rcolumn(fits, thdu, cols[COL_IFLUX].icol, 1, DAT_DBL, 1, NULL, 0, 1,
	       &sou->tot_flux) != 1)
      return 1;
/*
 * Get the source frequency offset for each IF.
 */
    for(i=0; i<ob->nif; i++) {
      if(rcolumn(fits, thdu, cols[COL_FREQOFF].icol, 1, DAT_DBL, 1, NULL, i, 1,
		 &ob->ifs[i].freq) != 1)
	return 1;
    };
  }
/*
 * No SU table? Attempt to ascertain source characteristics from the
 * information in the primary FITS header.
 */
  else {
/*
 * Get the source name.
 */
    if(phdu->object==NULL) {
      lprintf(stderr, "get_source: No source description in FITS file.\n");
      return 1;
    };
/*
 * Copy the source name into the Source descriptor.
 */
    stripcpy(sou->name, sizeof(sou->name), phdu->object, strlen(phdu->object));
/*
 * Get the RA and DEC and convert to radians.
 */
    sou->ra = dtor * get_axis(phdu, ax.rpos+1)->crval;
    sou->dec = dtor * get_axis(phdu, ax.dpos+1)->crval;
    sou->epoch = phdu->equinox;
/*
 * Assume FK4 if the epoch < 1984.0 (this is the convention introduced
 * in the 3rd wcs paper.
 */
    if(sou->epoch < 1984.0) {
      double ra2000,dec2000;   /* J2000 versions of the coordinates */
/*
 * Get the B1950.0 versions of the source coordinates.
 */
      double ra_B1950 = sou->ra;
      double dec_B1950 = sou->dec;
      if(fabs(sou->epoch-1950.0)>0.1)
	slaPreces("FK4", sou->epoch, 1950.0, &ra_B1950, &dec_B1950);
/*
 * Get the equivalent J2000 coordinates and use these to determine
 * the apparent RA/DEC.
 */
      slaFk45z(ra_B1950, dec_B1950, 1975.0, &ra2000, &dec2000);
      slaMap(ra2000, dec2000, 0.0, 0.0, 0.0, 0.0, 2000.0,
	     (date->ut + slaDat(date->utc_ref)) / daysec + date->utc_ref,
	     &sou->app_ra, &sou->app_dec);
/*
 * Assume FK5 if the coordinate system epoch is >= 1984.0.
 */
    } else {
      slaMap(sou->ra, sou->dec, 0.0, 0.0, 0.0, 0.0, sou->epoch,
	     (date->ut + slaDat(date->utc_ref)) / daysec + date->utc_ref,
	     &sou->app_ra, &sou->app_dec);
    };
/*
 * There are no further source details available.
 */
    sou->tot_flux = 0.0;
  };
/*
 * Attempt to find the observing center.
 */
  new_hline(fits->hdu, 0);   /* Rewind the primary header */
  while(next_key(fits, fits->hdu, misc_keys,
		 sizeof(misc_keys)/sizeof(Fitkey), EOH_SEEK, &key) == 0) {
    switch(key.keyid) {
    case OBSRA:
      obsra = KEYDBL(key) * dtor;
      have_obsra = 1;
      break;
    case OBSDEC:
      obsdec = KEYDBL(key) * dtor;
      have_obsdec = 1;
      break;
    };
  };
/*
 * If a pointing center was read, record it.
 */
  if(have_obsra && have_obsdec) {
    if(set_obs_radec(ob, obsra, obsdec))
      return 1;
/*
 * If the pointing center wasn't given, substitute the source location,
 * but mark this as tentative by setting have_obs to false, so that it
 * doesn't get recorded in output files.
 */
  } else {
    sou->have_obs = 0;
    sou->obsra = sou->ra;
    sou->obsdec = sou->dec;
    sou->east = 0.0;
    sou->north = 0.0;
  };
/*
 * List the source names.
 */
  lprintf(stdout, "Found source: %s\n", sou->name);
  return 0;
}

/*.......................................................................
 * Count the number of sources in an observation.
 *
 * Input:
 *  fob    Fitob *  The FITS-Observation intermediary descriptor.
 * Output:
 *  return   int    The number of sources counted, or 0 on error.
 */
static int count_sources(Fitob *fob)
{
  Thdu *thdu;  /* The SU table descriptor */
  int nsource; /* The number of sources */
/*
 * Look for an SU table.
 */
  thdu = find_table(fob->fits, "AIPS SU", 0, NULL);
/*
 * If there is no source table then ensure that the primary HDU header
 * cites a source.
 */
  if(thdu==NULL) {
    if(((Phdu *)fob->fits->hdu)->object==NULL) {
      lprintf(stderr, "count_sources: No source description in FITS file.\n");
      nsource = 0;
    } else {
      nsource = 1;
    };
  }
/*
 * Count the number of entries in the SU table.
 */
  else {
    nsource = numrow(thdu);
  };
  return nsource;
}

/*.......................................................................
 * Count the number of frequency groups in the FITS file.
 *
 * Input:
 *  fob    Fitob *  The FITS-Observation intermediary descriptor.
 * Output:
 *  return   int    The number of frequency groups counted, or 0 on error.
 */
static int count_FQ_entries(Fitob *fob)
{
  Thdu *thdu;  /* The FQ table descriptor */
/*
 * Look for an FQ table.
 */
  thdu = find_table(fob->fits, "AIPS FQ", 0, NULL);
  return thdu==NULL ? 1 : numrow(thdu);
}

/*.......................................................................
 * Determine and record the number of IFs in the FITS file, in fob->nif.
 *
 * Input:
 *  fob    Fitob *  The FITS/Observation intermediary descriptor.
 * Output:
 *  return   int    The number of IFs, or -1 on error.
 */
static int count_IFs(Fitob *fob)
{
/*
 * If there is an IF axis, then its dimension is the number of IFs.
 */
  fob->nif = (ax.ipos < 0) ? 1 : ((Phdu *) fob->fits->hdu)->dims[ax.ipos];
/*
 * Sanity check the number of IFs.
 */
  if(fob->nif <= 0) {
    lprintf(stderr, "count_IFs: Illegal IF axis dimension: %d\n", fob->nif);
    return -1;
  };
  return fob->nif;
}

/*.......................................................................
 * Read the first group and use its date to fill as the reference date
 * to use to fill fob->date.
 *
 * Input/Output:
 *  fob    Fitob *   The FITS/Observation structure being constructed.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int get_date(Fitob *fob)
{
  Parval *pval=NULL; /* Random group parameters */
  long gcount;       /* The number of groups in the file */
  int yy,mm,dd;      /* Gregorian Year, month and day of observation */
  int ierr=0;        /* Error return status from slalib routines */
  double day1;       /* UTC MJD of the first day of the observation */
  double year1;      /* UTC MJD of the start of the year of observation */
  double datutc;     /* Data time - UTC */
  double ddum;
/*
 * Search for the first visibility that has a valid date.
 * Note that on rare occasions UV data sets have been encountered that
 * contain one or more illegal dates of JD=0. When the UV file is
 * then sorted into TB order, these become the first visibilities
 * in the file. When it comes to binning and reading the UV data,
 * these visibilities will be ommitted.
 */
  gcount = fob->fits->hdu->gcount;
  for(fob->start_group=0; fob->start_group<gcount; fob->start_group++) {
/*
 * Read the next group.
 */
    pval = read_pars(fob, fob->start_group);
    if(!pval)
      return 1;
/*
 * Slalib slaDjcl requires mjd's of > -2395521.0 (1 March 4701BC).
 */
    if(pval->date > -2395521.0)
      break;
  };
/*
 * Have we reached the end of the data without finding a usable date?
 */
  if(!pval || fob->start_group >=gcount) {
    lprintf(stderr,
	    "get_date: There are no visibilities with valid dates.\n");
    return 1;
  };
/*
 * Warn about skipped integrations.
 */
  if(fob->start_group > 0L) {
    lprintf(stderr,
     "get_date: Skipped the first %ld visibilities. They had corrupt dates.\n",
	 fob->start_group);
  };
/*
 * Before we can interpret the start date, we need to know what time system
 * it belongs to. This potentially changes from one sub-array to the next
 * so it is recorded separately in each binary antenna table. In ASCII
 * tables and in binary tables that don't specify a time system, IAT is
 * the default.
 */
  if(get_subarray_time_systems(fob, slaDat(pval->date)))
    return 1;
/*
 * Use the antenna table that is accosiated with the integration
 * to deduce its time offset from UTC.
 */
  if(pval->isub<0 || pval->isub>fob->maxan) {
    lprintf(stderr, "get_date: Missing AN table, version: %d\n", pval->isub+1);
    return 1;
  };
  datutc = fob->antab[pval->isub].datutc;
/*
 * Get the date of the first day of the observation (UTC).
 */
  day1 = floor(pval->date - datutc/daysec);
/*
 * Get the year,month and day corresponding to day1.
 */
  slaDjcl(day1, &yy, &mm, &dd, &ddum, &ierr);
  if(ierr) {
    lprintf(stderr, "The first visibility has an unbelievable date (MJD=%g).\n",
	    day1);
    return 1;
  };
/*
 * Determine the UTC MJD corresponding to the start of the year of observation.
 */
  slaCldj(yy, 1, 1, &year1, &ierr);
  if(ierr) {
    lprintf(stderr, "Error translating the date of the first visibility.\n");
    return ierr;
  };
/*
 * Record the year.
 */
  fob->date.year = yy;
/*
 * Calculate and record the apparent GMST on day1.
 */
  fob->date.app_st = slaGmst(day1) + slaEqeqx(day1);
/*
 * Record year1 as the UTC origin from which to compute the number of days into
 * the year.
 */
  fob->date.utc_ref = year1;
/*
 * Calculate and record the number of seconds into the year of day1.
 */
  fob->date.ut = (day1 - fob->date.utc_ref) * daysec;
/*
 * Record the integration time that was determined in a prior call to
 * probe_times().
 */
  fob->date.cav_tim = fob->date.iav_tim = fob->binwid;
  return ierr;
}

/*.......................................................................
 * Read an AIPS AN table and record station details in the requested
 * sub-array descriptor.
 *
 * Input:
 *  ob  Observation * The descriptor of the observation being read.
 *  fits       Fits * The descriptor of the FITS file.
 *  an        Antab * AN table descriptor.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int get_stations(Observation *ob, Fits *fits, Antab *an)
{
  static const double vlalon = 1.878283678; /* VLA longitude (radians) */
  Subarray *sub;     /* Subarray descriptor */
  Thdu *thdu;        /* FITS descriptor of AN table */
  Phdu *phdu;        /* The primary header descriptor */
  Station *tel;      /* Pointer to station descriptor array */
  int i;
/*
 * Check args.
 */
  if(an==NULL || an->thdu==NULL) {
    lprintf(stderr, "get_stations: NULL antenna table descriptor.\n");
    return 1;
  };
/*
 * Get the subarray and FITS table descriptors for this antenna table.
 */
  sub = an->sub;
  thdu = an->thdu;
/*
 * Initialize the names in the station descriptor so that we can tell
 * if some don't get intitialized. Use the antwt member as an indication
 * of whether an entry for the antenna existed in the antenna table.
 */
  tel = sub->tel;
  for(i=0; i<sub->nstat; i++,tel++) {
    tel->name[0] = '\0';
    tel->type = GROUND;
    tel->antwt = 0;
  };
/*
 * The keywords in an ASCII 'AIPS AN' table are totally different to
 * those in a binary 'AIPS AN' table - so the two must be treated
 * differently - arggh.
 *
 * Binary or ascii 'AIPS AN' table?
 */
  switch(thdu->type) {
  case F_TABLE:
    if(getanasc(ob, fits, an))
      return 1;
    break;
  case F_BINTAB:
    if(getanbin(ob, fits, an))
      return 1;
    break;
  default:
    lprintf(stderr, "get_stations: AN table has unknown table type.\n");
    return 1;
  };
/*
 * Record the time offset of this sub-array.
 */
  sub->datutc = an->datutc;
/*
 * There are three known coordinate systems in use presently. The VLA
 * uses a VLA centered coordinate system, whereas VLBI and ACTA
 * use an Earth centered coordinate system. In addition the ACTA and
 * the VLA have the opposite signed Y axis. Make the required transformations
 * to the VLBI standard after deducing the array type from the PRIMARY
 * header TELESCOP keyword.
 */
  phdu = (Phdu *) fits->hdu;
  if(phdu->telescop) {
    if(strcmp(phdu->telescop, "VLA")==0) {
      for(i=0,tel = sub->tel; i<sub->nstat; i++,tel++) {
	if(tel->type == GROUND) {
	  tel->geo.gnd.x = tel->geo.gnd.x * cos(vlalon) +
	                   tel->geo.gnd.y * sin(vlalon);
	  tel->geo.gnd.y = tel->geo.gnd.x * sin(vlalon) -
	                   tel->geo.gnd.y * cos(vlalon);
	};
      };
    } else if(strcmp(phdu->telescop, "ACTA")==0) {
      for(i=0, tel = sub->tel; i<sub->nstat; i++,tel++) {
	if(tel->type == GROUND)
	  tel->geo.gnd.y = -tel->geo.gnd.y;
      };
    };
  };
/*
 * Binary AN tables record the telescope positions wrt a given
 * array center, which records the position wrt the center of the Earth in
 * meters. We require the absolute positions, so add in the array offset
 * here.
 */
  if(sub->binan) {
    Binan *binan = sub->binan;
    for(i=0, tel = sub->tel; i<sub->nstat; i++,tel++) {
      if(tel->type == GROUND) {
	tel->geo.gnd.x = binan->arrayx + tel->geo.gnd.x;
	tel->geo.gnd.y = binan->arrayy + tel->geo.gnd.y;
	tel->geo.gnd.z = binan->arrayz + tel->geo.gnd.z;
      };
    };
  };
/*
 * Check that all stations were initialized.
 */
  tel = sub->tel;
  for(i=0; i<sub->nstat; i++,tel++) {
    if(tel->antwt==0) {
      lprintf(stderr, "get_stations: Missing station %d in AN table %d.\n",
	      i+1, thdu->extver);
      return 1;
    };
/*
 * If the antenna has no name, generate a fake name for it.
 */
    if(tel->name[0] == '\0')
      sprintf(tel->name, "ANT%d", i+1);
  };
  return 0;
}

/*.......................................................................
 * Read a binary AIPS AN table and record station details in sub->tel[].
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation being read.
 *  fits       Fits *  The descriptor of the FITS file.
 *  an        Antab *  The intermediary AN descriptor.
 * Output:
 *  return      int     0 - OK.
 *                      1 - Error.
 */
static int getanbin(Observation *ob, Fits *fits, Antab *an)
{
  Thdu *thdu;        /* FITS AN table descriptor. */
  Subarray *sub;     /* The subarray to be processed. */
  Fitkey key;        /* Descriptor of FITS keyword/value pair */
  double iatutc=0.0; /* Archaic alternate for DATUTC */
  Binan *binan;      /* Container for FITS style sub-array info */
  int nchar;         /* Number of characters read from char column */
  int icol;          /* Column number - 1-relative */
  int irow;          /* Row number - 1-relative */
  int iant;          /* Telescope index 0-relative */
  int nopcal;        /* Number of polarization cal parameters in the table */
  int numorb;        /* Number of orbital cal parameters in the table */
  int i;
/*
 * Create a table of optional keywords.
 */
  enum {ARRAYX, ARRAYY, ARRAYZ, GSTIA0, DEGPDY, AN_FREQ, RDATE,	POLARX, POLARY,
	UT1UTC, DATUTC, IATUTC, TIMSYS, ARRNAM, NUMORB, POLTYPE};
  static Fitkey ankeys[]={
    {"ARRAYX", 0, ARRAYX,   DAT_DBL, NULL, NULL},
    {"ARRAYY", 0, ARRAYY,   DAT_DBL, NULL, NULL},
    {"ARRAYZ", 0, ARRAYZ,   DAT_DBL, NULL, NULL},
    {"FREQ",   0, AN_FREQ,  DAT_DBL, NULL, NULL},
    {"GSTIA0", 0, GSTIA0,   DAT_DBL, NULL, NULL},
    {"DEGPDY", 0, DEGPDY,   DAT_DBL, NULL, NULL},
    {"RDATE",  0, RDATE,    DAT_STR, NULL, NULL},
    {"POLARX", 0, POLARX,   DAT_DBL, NULL, NULL},
    {"POLARY", 0, POLARY,   DAT_DBL, NULL, NULL},
    {"UT1UTC", 0, UT1UTC,   DAT_DBL, NULL, NULL},
    {"DATUTC", 0, DATUTC,   DAT_DBL, NULL, NULL},
    {"IATUTC", 0, IATUTC,   DAT_DBL, NULL, NULL},
    {"TIMSYS", 0, TIMSYS,   DAT_STR, NULL, NULL},
    {"ARRNAM", 0, ARRNAM,   DAT_STR, NULL, NULL},
    {"POLTYPE",0, POLTYPE,  DAT_STR, NULL, NULL},
  };
/*
 * Define a structure used to record field assignments.
 */
  typedef struct {
    char *name;    /* Name of field */
    int icol;      /* Number of field column in row (0-relative) */
    int need;      /* True if this column is mandatory */
  } Anfield;
/*
 * Define BINTAB-table-style 'AIPS AN' field assignments.
 */
  enum {ANNAME, NOSTA, MNTSTA, STABXYZ, ORBPARM, STAXOF,
        POLTYA, POLAA, POLCALA,	POLTYB, POLAB, POLCALB};  /* fields[] indexes */
  static Anfield fields[]={
    {"ANNAME",  0, 1},
    {"NOSTA",   0, 1},
    {"MNTSTA",  0, 1},
    {"STABXYZ", 0, 1},
    {"ORBPARM", 0, 1},
    {"STAXOF",  0, 0},
    {"POLTYA",  0, 0},
    {"POLAA",   0, 0},
    {"POLCALA", 0, 0},
    {"POLTYB",  0, 0},
    {"POLAB",   0, 0},
    {"POLCALB", 0, 0},
  };
  static const int nfield=sizeof(fields)/sizeof(Anfield);
  Anfield *field;   /* Pointer into fields[] */
/*
 * Get the AN table descriptor and output sub-array descriptor.
 */
  thdu = an->thdu;
  sub = an->sub;
/*
 * Get the value of the keyword that specifies the number of orbital
 * parameters recorded for each antenna.
 */
  if(get_key(fits, (Hdu *)thdu, "NUMORB", DAT_INT, LOOP_SEEK, &key)==KEY_FOUND)
    numorb = KEYINT(key);
  else
    numorb = 0;
/*
 * Get the value of the keyword that specifies the number of polarization
 * cal parameters recorded for each feed of each antenna.
 */
  if(get_key(fits, (Hdu *)thdu, "NOPCAL", DAT_INT, LOOP_SEEK, &key)==KEY_FOUND)
    nopcal = KEYINT(key);
  else
    nopcal = 0;
/*
 * Allocate a Binan descriptor in the current sub-array to record the
 * contents of the table. They are preserved their-in for copying to
 * output FITS files later.
 */
  binan = new_Binan(sub, sub->nstat, nopcal, numorb);
  if(binan==NULL)
    return 1;
/*
 * Read the rest of the expected binary table AN keyword/value pairs.
 */
  new_hline((Hdu *) thdu, 0);   /* Rewind header */
  while(next_key(fits, (Hdu *) thdu, ankeys,
		 sizeof(ankeys)/sizeof(Fitkey), EOH_SEEK, &key) == 0) {
    switch(key.keyid) {
    case ARRAYX:
      binan->arrayx = KEYDBL(key);
      break;
    case ARRAYY:
      binan->arrayy = KEYDBL(key);
      break;
    case ARRAYZ:
      binan->arrayz = KEYDBL(key);
      break;
    case GSTIA0:
      binan->gstia0 = KEYDBL(key);
      break;
    case DEGPDY:
      binan->degpdy = KEYDBL(key);
      break;
    case AN_FREQ:
      binan->freq = KEYDBL(key);
      break;
    case RDATE:
      stripcpy(binan->rdate, sizeof(binan->rdate), KEYSTR(key),
	       strlen(KEYSTR(key)));
      break;
    case POLARX:
      binan->polarx = KEYDBL(key);
      break;
    case POLARY:
      binan->polary = KEYDBL(key);
      break;
    case UT1UTC:
      binan->ut1utc = KEYDBL(key);
      break;
    case DATUTC:
      binan->datutc = KEYDBL(key);
      break;
    case IATUTC:
      iatutc = KEYDBL(key);
      break;
    case TIMSYS:
      stripcpy(binan->timsys, sizeof(binan->timsys), KEYSTR(key),
	       strlen(KEYSTR(key)));
      break;
    case ARRNAM:
      stripcpy(binan->arrnam, sizeof(binan->arrnam), KEYSTR(key),
	       strlen(KEYSTR(key)));
      break;
    case POLTYPE:
      stripcpy(binan->poltype, sizeof(binan->poltype), KEYSTR(key),
	       strlen(KEYSTR(key)));
      break;
    };
  };
/*
 * Read the reference antenna number and associated R-L phase differences
 * if present.
 */
  if(rd_p_refant(fits, (Hdu *)thdu, sub))
    return 1;
/*
 * The default for TIMSYS is "IAT" according to Going AIPS.
 */
  if(string_is_empty(binan->timsys))
    strcpy(binan->timsys, "IAT");
/*
 * The correct value for datutc was previously determined by
 * get_subarray_time_systems() and recorded in an->datutc.
 */
  binan->datutc = an->datutc;
/*
 * Search for each of the fields named in fields[].
 */
  for(field=fields,i=0; i<nfield; i++,field++) {
    field->icol = find_column(thdu, field->name, 0);
    if(field->icol==0 && field->need) {
      lprintf(stderr, "getanbin: Missing %s field in AN table.\n", field->name);
      return 1;
    };
  };
/*
 * Read entries for each antenna in the table.
 */
  for(irow=1; irow<=an->nrow; irow++) {
    Anmap *smap;            /* Station map descriptor */
/*
 * Read the 1-relative station number.
 */
    if(rcolumn(fits, thdu, fields[NOSTA].icol, irow, DAT_INT, 1, NULL, 0, 1,
	       &iant)!=1)
      return -1;
/*
 * Get the station map descriptor corresponding to antenna number 'iant'.
 */
    smap = &an->smap[an->antrow[iant]];
/*
 * Do we want this station?
 */
    if(smap->used) {
/*
 * Get the output telescope descriptor, and the associated
 * member of the descriptor used to record binary table specific parameters.
 */
      Station *tel = &sub->tel[smap->slot];
      Bintel *bt = &binan->bt[smap->slot];
/*
 * Read the telescope name.
 */
      nchar = rcolumn(fits, thdu, fields[ANNAME].icol, irow, DAT_CHR, 1, NULL,
		      0, (int) sizeof(bt->anname), &bt->anname);
      if(nchar < 1)
	return -1;
      stripstr(bt->anname, nchar);
/*
 * Copy the station name into the telescope descriptor.
 */
      stripcpy(tel->name, sizeof(tel->name), bt->anname, strlen(bt->anname));
/*
 * Record the antenna number.
 */
      tel->antno = iant;
      bt->nosta = iant;
/*
 * Get the mount type.
 */
      if(rcolumn(fits, thdu, fields[MNTSTA].icol, irow, DAT_INT, 1, NULL,
		 0, 1, &bt->mntsta)!=1)
	return 1;
/*
 * Get the array coordinates.
 */
      if(rcolumn(fits, thdu, fields[STABXYZ].icol, irow, DAT_DBL, 1, NULL,
		 0, 3, bt->stabxyz) != 3L)
	return 1;
/*
 * Get orbital parameters if present.
 */
      if(binan->numorb > 0) {
	if(rcolumn(fits, thdu, fields[ORBPARM].icol, irow, DAT_DBL, 1,
		   NULL, 0, binan->numorb, bt->orbparm) != binan->numorb)
	  return 1;
      };
/*
 * Read the X,Y,Z or orbital parameter arrays and record
 * in the telescope descriptors.
 */
      switch(bt->mntsta) {
      default:
	if(!(strcmp(binan->arrnam, "CBI")==0 && bt->mntsta==4)) {
	  lprintf(stderr, "Warning: Unknown AN-table MNTSTA value (%d).\n",
		  bt->mntsta);
	  lprintf(stderr, "         Will assume that it is ground-based.\n");
	};
/* Note fallthrough to case 0 */
      case 0:  /* altaz */
      case 1:  /* equatorial */
	tel->type = GROUND;
	tel->geo.gnd.x = bt->stabxyz[0];
	tel->geo.gnd.y = bt->stabxyz[1];
	tel->geo.gnd.z = bt->stabxyz[2];
	break;
      case 2:
	tel->type = ORBIT;
	if(binan->numorb >= 6) {
	  tel->geo.orb.semi_major = bt->orbparm[0];
	  tel->geo.orb.eccentricity = bt->orbparm[1];
	  tel->geo.orb.inclination = bt->orbparm[2];
	  tel->geo.orb.ra_ascending = bt->orbparm[3];
	  tel->geo.orb.arg_perigee = bt->orbparm[4];
	  tel->geo.orb.mean_anomoly = bt->orbparm[5];
	};
      };
/*
 * Get the rest of the parameters, to be recorded in bt only.
 */
/*
 * Axis offset.
 */
      icol = fields[STAXOF].icol;
      if(icol>0 && rcolumn(fits, thdu, icol, irow, DAT_DBL, 1, NULL, 0,
			   1, &bt->staxof)!=1)
	return 1;
/*
 * Feed A polarization type.
 */
      icol = fields[POLTYA].icol;
      if(icol>0 && rcolumn(fits, thdu, icol, irow, DAT_CHR, 1, NULL, 0,
			   1, &bt->poltya)!=1)
	return 1;
/*
 * Feed A position angle.
 */
      icol = fields[POLAA].icol;
      if(icol>0 && rcolumn(fits, thdu, icol, irow, DAT_DBL, 1, NULL, 0,
			   1, &bt->polaa)!=1)
	return 1;
/*
 * Feed A polarization cal parameters.
 */
      icol = fields[POLCALA].icol;
      if(icol>0 && rcolumn(fits, thdu, icol, irow, DAT_DBL, 1, NULL, 0,
			   binan->nopcal, bt->polcala)!=binan->nopcal)
	return 1;
/*
 * Feed B polarization type.
 */
      icol = fields[POLTYB].icol;
      if(icol>0 && rcolumn(fits, thdu, icol, irow, DAT_CHR, 1, NULL, 0,
			   1, &bt->poltyb)!=1)
	return 1;
/*
 * Feed B position angle.
 */
      icol = fields[POLAB].icol;
      if(icol>0 && rcolumn(fits, thdu, icol, irow, DAT_DBL, 1, NULL, 0,
			   1, &bt->polab)!=1)
	return 1;
/*
 * Feed B polarization cal parameters.
 */
      icol = fields[POLCALB].icol;
      if(icol>0 && rcolumn(fits, thdu, icol, irow, DAT_DBL, 1, NULL, 0,
			   binan->nopcal, bt->polcalb)!=binan->nopcal)
	return 1;
/*
 * Mark the antenna as now defined.
 */
      tel->antwt = 1.0;
    };
  };
  return 0;
}

/*.......................................................................
 * Read an ASCII AIPS AN table and record station details in sub->tel[].
 *
 * Input:
 *  ob  Observation * The descriptor of the observation being read.
 *  fits       Fits * The descriptor of the FITS file.
 *  an        Antab *  The intermediary AN descriptor.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int getanasc(Observation *ob, Fits *fits, Antab *an)
{
  Thdu *thdu;      /* FITS AN table descriptor. */
  Subarray *sub;   /* The subarray to be processed. */
  int nchar;       /* Number of characters read from char column */
  int irow;        /* Row number - 1-relative */
  int iant;        /* Telescope index 0-relative */
  int i;
/*
 * Define a structure used to record field assignments.
 */
  typedef struct {
    char *name;    /* Name of field */
    int icol;      /* Number of field column in row (0-relative) */
    int need;      /* True if this axis is mandatory */
  } Anfield;
/*
 * Define ASCII-table-style 'AIPS AN' field assignments.
 */
  enum {ANT_NO, STATION, LX, LY, LZ};  /* Used to index 'fields[]' */
  static Anfield fields[]={
    {"ANT NO.", 0, 1},
    {"STATION", 0, 1},
    {"LX",      0, 1},
    {"LY",      0, 1},
    {"LZ",      0, 1},
  };
  static const int nfield=sizeof(fields)/sizeof(Anfield);
  Anfield *field;   /* Pointer into fields[] */
/*
 * Get the AN table descriptor and output sub-array descriptor.
 */
  thdu = an->thdu;
  sub = an->sub;
/*
 * Read the reference antenna number and associated R-L phase differences
 * if present.
 */
  if(rd_p_refant(fits, (Hdu *)thdu, sub))
    return 1;
/*
 * Search for each of the fields named in fields[].
 */
  for(field=fields,i=0; i<nfield; i++,field++) {
    field->icol = find_column(thdu, field->name, 0);
    if(field->icol==0 && field->need) {
      lprintf(stderr, "getanasc: Missing %s field in AN table.\n", field->name);
      return 1;
    };
  };
/*
 * Read entries for each antenna in the table.
 */
  for(irow=1; irow<=an->nrow; irow++) {
    Anmap *smap;            /* Station map descriptor */
/*
 * Read the 1-relative station number.
 */
    if(rcolumn(fits, thdu, fields[ANT_NO].icol, irow, DAT_INT, 1, NULL, 0, 1,
	       &iant)!=1)
      return -1;
/*
 * Get the station map descriptor corresponding to antenna number 'iant'.
 */
    smap = &an->smap[an->antrow[iant]];
/*
 * Do we want this station?
 */
    if(smap->used) {
/*
 * Get the output telescope descriptor.
 */
      Station *tel = &sub->tel[smap->slot];
/*
 * Read the telescope name.
 */
      nchar = rcolumn(fits, thdu, fields[STATION].icol, irow, DAT_CHR, 1, NULL,
		      0, (int) sizeof(tel->name), &tel->name);
      if(nchar < 1)
	return -1;
      stripstr(tel->name, nchar);  /* Remove trailing spaces */
/*
 * Record the antenna number.
 */
      tel->antno = iant;
/*
 * Read the X,Y, and Z fields and record in the telescope descriptor.
 */
      if(rcolumn(fits, thdu, fields[LX].icol, irow, DAT_DBL, 1, NULL, 0,
		 1, &tel->geo.gnd.x)!=1L)
	return 1;
      if(rcolumn(fits, thdu, fields[LY].icol, irow, DAT_DBL, 1, NULL, 0,
		 1, &tel->geo.gnd.y)!=1L)
	return 1;
      if(rcolumn(fits, thdu, fields[LZ].icol, irow, DAT_DBL, 1, NULL, 0,
		 1, &tel->geo.gnd.z)!=1L)
	return 1;
/*
 * Mark the antenna as now defined.
 */
      tel->antwt = 1.0;
    };
  };
  return 0;
}

/*.......................................................................
 * Determine baseline parameters and record them in sub->base[].
 *
 * This function requires that the following members have been initialized:
 * In fob->an[*]:
 *  bmap[],nbase,nstat
 * In ob:
 *  sub[isub].tel[], source, app_st
 * It also requires that the baseline array ob->sub[isub].base[] has been
 * allocated for ob->sub[isub].nbase baselines.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation being read.
 *  fob       Fitob *  The FITS/Observation intermediary descriptor.
 *  an        Antab *  AN table descriptor.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int get_baselines(Observation *ob, Fitob *fob, Antab *an)
{
  Subarray *sub;     /* Subarray descriptor */
  Baseline *b;       /* Pointer to baseline descriptor in ob->base[] */
  Station *tel_a, *tel_b; /* Pointers to telescope descriptors */
  double bx, by, bz; /* Baseline components */
  int ia,ib;         /* Telescope numbers */
  Anmap *bmap;       /* Pointer to entry in an->bmap[] */
  Anmap *smap;       /* Pointer to an->smap[] */
/*
 * Get the Subarray descriptor for this AN table.
 */
  sub = an->sub;
/*
 * Loop through all potential baselines and record each used baseline
 * (as determined from the associated flag in fob->bmap[]) in the next
 * available baseline descriptor of sub->base[].
 */
  smap = an->smap;
  bmap = an->bmap;
  for(ia=0; ia<an->nrow; ia++) {
    for(ib=ia+1; ib<an->nrow; ib++,bmap++) {
/*
 * Is this baseline required?
 */
      if(bmap->used) {
/*
 * Get the baseline descriptor mapped to the baseline.
 */
	b = &sub->base[bmap->slot];
/*
 * Map from AN table telescope indexes to sub->tel[] indexes.
 */
	b->tel_a = smap[ia].slot;
	b->tel_b = smap[ib].slot;
/*
 * Get the respective telescope descriptors.
 */
	tel_a = &sub->tel[b->tel_a];
	tel_b = &sub->tel[b->tel_b];
/*
 * We can only support this for ground-based telescopes.
 */
	if(tel_a->type==GROUND && tel_b->type==GROUND) {
/*
 * Get relative distances between telescopes.
 */
	  bx = (tel_a->geo.gnd.x - tel_b->geo.gnd.x);
	  by = (tel_a->geo.gnd.y - tel_b->geo.gnd.y);
	  bz = (tel_a->geo.gnd.z - tel_b->geo.gnd.z);
/*
 * Compute and install the baseline, hour-angle, distance-in-meters,
 * and Z-distance-in-meters.
 */
	  b->boff  = ob->date.app_st - ob->source.app_ra - pi/2.0 - 
	    (by==0.0 && bx==0.0) ? 0.0 : atan2(by,bx);
	  b->bxy   = sqrt(bx*bx + by*by);
	  b->bz    = bz;
	};
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Determine IF-specific info such as IF frequencies and record them
 * in ob->ifs[]. The following functions need to be called before
 * this function:
 *  Obs_alloc:     Creates ob.
 *  get_source:    Adds source freq. offsets in ob->ifs[*].freq .
 *
 * Input:
 *  ob Observation *  The descriptor of the Observation being read.
 *  fits      Fits *  The descriptor of the FITS file.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int get_IF_freq(Observation *ob, Fits *fits)
{
  Phdu *phdu;        /* Pointer to primary HDU descriptor */
  Thdu *thdu;        /* The FQ-table descriptor */
  Imaxis *axis;      /* The FITS descriptor of the FREQ axis */
  If *ifptr;         /* POinter to IF descriptor being initialized */
  int i;
/*
 * Get a pointer to the primary HDU descriptor.
 */
  phdu = (Phdu *) fits->hdu;
/*
 * Get the descriptor for the FREQ primary HDU axis.
 */
  axis = get_axis(phdu, ax.fpos+1);
  if(axis==NULL)
    return 1;
/*
 * Look for an FQ table.
 */
  thdu = find_table(fits, "AIPS FQ", 0, NULL);
/*
 * No FQ table?
 */
  if(thdu==NULL) {
/*
 * This is an error if there is an IF axis.
 */
    if(ax.ipos >= 0 || ob->nif > 1) {
      lprintf(stderr, "get_IF_freq: Unable to locate AIPS FQ table.\n");
      return 1;
    };
/*
 * Get the single IF descriptor to be initialized.
 */
    ifptr = &ob->ifs[0];
/*
 * No FQ table and no IF axis. The required base frequency is that of the
 * first pixel on the FREQ axis - add this to the already assigned source
 * frequency offset.
 */
    ifptr->freq += axis->crval + axis->cdelt * (1.0 - axis->crpix);
    ifptr->df = axis->cdelt;
    ifptr->bw = fabs(ifptr->df) * phdu->dims[ax.fpos];
    ifptr->coff = 0;
  }
/*
 * An FQ table was found.
 */
  else {
/*
 * Define a structure used to record field assignments.
 */
    typedef struct {
      char *name;    /* Name of field */
      int icol;      /* Number of field column in row (0-relative) */
      int need;      /* True if this axis is mandatory */
    } FQfield;
/*
 * Define 'AIPS FQ' field assignments.
 */
    enum {IF_FREQ, CH_WIDTH, TOTAL_BW, SIDEBAND}; /* Used to index 'fields[]'*/
    static FQfield fields[]={
      {"IF FREQ",         0, 1},
      {"CH WIDTH",        0, 1},
      {"TOTAL BANDWIDTH", 0, 1},
      {"SIDEBAND",        0, 1},
    };
    static const int nfield=sizeof(fields)/sizeof(FQfield);
    FQfield *field;   /* Pointer into fields[] */
/*
 * Search for each of the fields named in fields[].
 */
    for(field=fields,i=0; i<nfield; i++,field++) {
      field->icol = find_column(thdu, field->name, 0);
      if(field->icol==0 && field->need) {
	lprintf(stderr, "get_IF_freq: Missing %s field in FQ table.\n",
		field->name);
	return 1;
      };
    };
/*
 * If there is more than one row in the table, signal an error.
 * We can't handle changes in IF frequencies definitions.
 */
    if(numrow(thdu) != 1) {
      lprintf(stderr, "get_IF_freq: The FQ table has more than 1 FREQID - this can not be handled.\n");
      return 1;
    };
/*
 * Read the required members of the table row for each IF.
 */
    ifptr = &ob->ifs[0];
    for(i=0; i<ob->nif; i++,ifptr++) {
      double dtmp;
      int itmp;
/*
 * Get the IF frequency offset and add it to the already assigned source
 * frequency offset for this IF.
 */
      if(rcolumn(fits, thdu, fields[IF_FREQ].icol, 1, DAT_DBL, 1, NULL, i,
		 1, &dtmp)!=1L)
	return 1;
      ifptr->freq += dtmp;
/*
 * Get the spectral-line channel width in this IF.
 */
      if(rcolumn(fits, thdu, fields[CH_WIDTH].icol, 1, DAT_DBL, 1, NULL, i,
		 1, &dtmp)!=1L)
	return 1;
      ifptr->df = dtmp;
/*
 * Get the total bandwidth of the IF.
 */
      if(rcolumn(fits, thdu, fields[TOTAL_BW].icol, 1, DAT_DBL, 1, NULL, i,
		 1, &dtmp)!=1L)
	return 1;
      ifptr->bw = dtmp;
/*
 * Get the sideband type recorded in this IF.
 */
      if(rcolumn(fits, thdu, fields[SIDEBAND].icol, 1, DAT_INT, 1, NULL, i,
		 1, &itmp)!=1L)
	return 1;
      if(itmp < 0)
	ifptr->df = -fabs(ifptr->df);
/*
 * Now add in the base frequency for channel 1.
 */
      ifptr->freq += axis->crval + ifptr->df * (1.0 - axis->crpix);
/*
 * Also record the channel offset for the new IF.
 */
      ifptr->coff = i * ob->nchan;
    };
  };
/*
 * Report number of IFs and spectral-line channels.
 */
  lprintf(stdout, "\nThere %s %d IF%s, and a total of %d channel%s:\n",
	  ob->nif>1?"are":"is", ob->nif, ob->nif>1?"s":"",
	  ob->nctotal, ob->nctotal>1?"s":"");
/*
 * Report IF characteristics.
 */
  lprintf(stdout, "\n %s\n %s\n %s\n",
    "IF  Channel    Frequency  Freq offset  Number of   Overall IF",
    "     origin    at origin  per channel   channels    bandwidth",
    "------------------------------------------------------------- (Hz)");
  for(i=0; i<ob->nif; i++) {
    If *ifptr = &ob->ifs[i];
    lprintf(stdout, " %2.2d  %7d %12g %12g    %7d %12g\n", i+1, ifptr->coff+1,
	    ifptr->freq, ifptr->df, ob->nchan, ifptr->bw);
  };
  return 0;
}

/*.......................................................................
 * Determine and record the number of stokes parameters or polarizations
 * in the FITS file, in fob->npol.
 *
 * Input:
 *  fob    Fitob *  The FITS/Observation intermediary descriptor.
 * Output:
 *  return   int    The number of polarizations, or -1 on error.
 */
static int count_stokes(Fitob *fob)
{
/*
 * If there is a STOKES axis, then its dimension is the number
 * of polarizations or stokes parameters recorded.
 */
  fob->npol = (ax.spos < 0) ? 1 : ((Phdu *) fob->fits->hdu)->dims[ax.spos];
/*
 * Sanity check the number of polarizations.
 */
  if(fob->npol <= 0) {
    lprintf(stderr, "count_stokes: Illegal STOKES axis dimension: %d\n",
	    fob->npol);
    return -1;
  };
  return fob->npol;
}

/*.......................................................................
 * Determine the types of polarizations or stokes parameters that were
 * observed and these types in ob->pols[].
 * The following functions need to be called before this function:
 *  count_stokes:  Determines the number of polarizations recorded (npol).
 *  Obs_alloc:     Creates ob,ob->npol and ob->pols[1..npol].
 *
 * Input:
 *  ob Observation *  The descriptor of the Observation being read.
 *  fits      Fits *  The descriptor of the FITS file.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int get_stokes(Observation *ob, Fits *fits)
{
  int i;
/*
 * If there is no polarization axis then assign I polarization.
 */
  if(ax.spos < 0) {
    ob->pols[0] = SI;
  }
/*
 * If there is a STOKES axis, then the polarizations are encoded with
 * the same numbers assigned in the Stokes enumeration, incrementing from
 * a reference value by an equal amount per stokes element on the stokes axis.
 */
  else {
    Imaxis *axis; /* Pointer to stokes axis descriptor */
    int ipol;     /* FITS polarization enumerator */
/*
 * Get a pointer to the stokes axis descriptor.
 */
    axis = get_axis((Phdu *)fits->hdu, ax.spos+1);
/*
 * Check and assign polarization types into ob->pols[].
 */
    for(i=0; i<ob->npol; i++) {
      ipol = axis->crval + ((i+1) - axis->crpix) * axis->cdelt;
      switch(ipol) {
      case SI: case SQ: case SU: case SV:
      case RR: case LL: case RL: case LR:
      case XX: case YY: case XY: case YX:
	ob->pols[i] = (Stokes) ipol;
	break;
      default:
	lprintf(stderr, "get_stokes: Unknown stokes type enumerator: %d\n",
		ipol);
	return 1;
	break;
      };
    };
  };
/*
 * Report the polarizations found.
 */
  lprintf(stdout, "\nPolarization(s):");
  for(i=0; i<ob->npol; i++)
    lprintf(stdout, " %s", Stokes_name(ob->pols[i]));
  lprintf(stdout, "\n");
  return 0;
}

/*.......................................................................
 * Read and record FITS history lines from the header of the primary HDU.
 * The lines are recorded in ob->his. Certain special keywords may
 * be encoded in history lines and these will be detected and recorded
 * where required.
 *
 * Input:
 *  ob Observation *  The descriptor of the Observation being read.
 *  fits      Fits *  The descriptor of the FITS file.
 * Input/Output:
 *  fob      Fitob *  If the AIPS WTSCAL parameter is detected, its
 *                    absolute value will be stored in ob->geom.wtscale
 *                    and its sign in fob->wtsign.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
static int get_history(Observation *ob, Fits *fits, Fitob *fob)
{
  Hdu *hdu = fits->hdu; /* Primary HDU descriptor */
  Fitkey key;           /* Keyword value pair returned by get_key() */
/*
 * Rewind to the first line of the primary HDU.
 */
  new_hline(hdu, 0);
/*
 * Read and record one history line at a time.
 */
  while(get_key(fits, hdu, "HISTORY", DAT_COM, EOH_SEEK, &key) == KEY_FOUND) {
/*
 * Filter out special AIPS HISTORY lines.
 */
    check_history(ob, KEYSTR(key), fob);
/*
 * Record the history line in the history.scr scratch file.
 */
    if(add_hist(ob, KEYSTR(key)))
      return 1;
  };
/*
 * Report the number of lines read.
 */
  lprintf(stdout, "\nRead %d lines of history.\n", ob->nhist);
  return 0;
}

/*.......................................................................
 * AIPS encodes a few of its keywords in FITS history lines.
 * Check for relevant lines here and extract values where relevant.
 *
 * Note that AIPS retains all history lines. This means that the same
 * special HISTORY lines often appear more than once. The last
 * such instance is the most recent and the one to be used. Thus make
 * sure that each time a recognised line is found, its value
 * over-rides any previous value of the same form.
 *
 * Recorded special values include:
 *
 * AIPS WTSCAL :   Absolute value recorded in ob->geom.wtscale,
 *                 sign recorded in fob->wtsign.
 *
 * Input:
 *  ob     Observation *  The descriptor of the observation.
 *  hline         char *  The text of the HISTORY line.
 *  fob          Fitob *  The FITS/Observation intermediary descriptor.
 *                        Some special values will be returned in this
 *                        descriptor - currently only fob->wtsign.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
static int check_history(Observation *ob, char *hline, Fitob *fob)
{
  int i;
/*
 * Skip white-space to the first word in the text of the history line.
 */
  while(isspace((int)*hline))
    hline++;
/*
 * Special AIPS history lines all start with the word AIPS.
 */
  if(strncmp(hline, "AIPS ", 5)==0) {
    int type;  /* The enumerated type of keyword found */
/*
 * List special AIPS history keywords.
 */
    typedef struct {
      char *name;
      enum {WTSCAL} type;
    } Aipskey;
    Aipskey aipskeys[] = {
      {"WTSCAL", WTSCAL},
    };
/*
 * Skip the "AIPS " prefix.
 */
    hline += 5;
/*
 * Skip white-space to the special KEYWORD name.
 */
    while(isspace((int)*hline))
      hline++;
/*
 * Check for recognised AIPS keywords.
 */
    type = -1;
    for(i=0; i<sizeof(aipskeys)/sizeof(aipskeys[0]); i++) {
      Aipskey *key = &aipskeys[i];
      int slen = strlen(key->name);
/*
 * A keyword should be followed by white-space or '='.
 */
      if(strncmp(hline, key->name, slen) == 0 &&
	 (isspace((int)hline[slen]) || hline[slen] == '=')) {
/*
 * Keyword recognised - record its enumerated identifier and
 * advance the hline pointer past the keyword.
 */
	type = (int) key->type;
	hline += slen;
	break;
      };
    };
/*
 * If the keyword was recognised prepare to read the associated value.
 */
    if(type >= 0) {
/*
 * Skip any white-space and '=' characters preceding the value.
 */
      while(isspace((int)*hline) || *hline=='=')
	hline++;
/*
 * The AIPS WTSCAL parameter contains a scale factor to be applied
 * to the weights.
 */
      switch(type) {
      case WTSCAL:
/*
 * Read and record the scale factor for later application.
 * This will overwrite the values from any previous instances of the
 * AIPS WTSCAL line.
 */
	{
	  double wtscale = atof(hline);
/*
 * Record the sign separately for direct application to the data.
 */
	  fob->wtsign = wtscale<0.0 ? -1:1;
/*
 * Record the absolute value of the scale factor. This will be applied
 * only as data are paged into memory.
 */
	  ob->geom.wtscale = fabs(wtscale);
	  if(ob->geom.wtscale == 0.0)
	    ob->geom.wtscale = 1;
	};
	break;
      default:
	break;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Read the data section of a given group from the FITS file.
 *
 * Input:
 *  fob    Fitob *  The FITS Observation intermediary descriptor.
 *  group   long    The 0-relative number of the group to read.
 * Output:
 *  pars  double *  Pointer to fob->data which contains fob->ndata
 *                  group parameters, or NULL on error.
 */
static double *get_data(Fitob *fob, long group)
{
  if(rimage(fob->fits, (Phdu *)fob->fits->hdu, group, 0L, fob->ndata, DAT_DBL,1,
	    NULL, fob->data) != fob->ndata)
    return NULL;
  else
    return fob->data;
}

/*.......................................................................
 * Given two FITS station indexes, return a pointer to the
 * element of the appropriate fob->antab[pval->isub].bmap[] for the
 * baseline of stations pval->ta, pval->tb.
 *
 * Input:
 *  fob         Fitob * The FITS/Observation intermediary descriptor.
 *  pval       Parval * Random-group parameters from read_pars().
 * Output:
 *  return      Anmap * The pointer to the associated baseline map, or
 *                      NULL if the baseline should be discarded.
 *                      Examples of baselines that are to be discarded,
 *                      are those of autocerrlation visibilities, where
 *                      ta==tb, and those in which one or both of the
 *                      stations doesn't exist in the respective
 *                      antenna table.
 */
static Anmap *loc_base(Fitob *fob, Parval *pval)
{
  static int ib;    /* Index of fob->bmap[] element */
  static int isub;  /* AN sub-array index */
  static int ta,tb; /* Antenna table antenna row-indexes */
  static int nstat; /* Number of stations in sub-array */
  static Antab *an; /* Antenna table indexing descriptor */
/*
 * Get the AN table extension number.
 */
  isub = pval->isub;
  if(isub<0 || isub>fob->maxan) {
    lprintf(stderr, "loc_base: Missing AN table, version: %d\n", isub+1);
    return NULL;
  };
/*
 * Get the Antab descriptor for the AN table.
 */
  an = &fob->antab[isub];
/*
 * No bmap array.
 */
  if(an->bmap == NULL) {
    lprintf(stderr, "loc_base: Missing AN table, version: %d\n", isub+1);
    return NULL;
  };
/*
 * Get the number of stations in this sub-array.
 */
  nstat = an->nrow;
/*
 * Check the telescope numbers.
 */
  if(pval->ta<0 || pval->ta>an->nsmax ||
     pval->tb<0 || pval->tb>an->nsmax)  {
    lprintf(stderr,"loc_base: Telescope index(es) %d,%d in data, but not in AN table.\n", pval->ta, pval->tb);
    return NULL;
  };
/*
 * Ignore autocorrelation visibilities.
 */
  if(pval->ta == pval->tb)
    return NULL;
/*
 * Convert from telescope number to AN table row index.
 */
  ta = an->antrow[pval->ta];
  tb = an->antrow[pval->tb];
/*
 * Swap ta and tb if ta > tb.
 */
  if(ta > tb) {int tmp=ta; ta=tb; tb=tmp;};
/*
 * Are the telescope indexes reasonable?
 */
  if(ta < 0 || tb >= an->nrow) {
    lprintf(stderr,"loc_base: Telescope index in data, but not in AN table.\n");
    return NULL;
  };
/*
 * The associated element in fob->bmap[] is defined as follows.
 * For increasing values of ta from 0 to nstat-1, ta baselines are
 * recorded for tb values from ta+1 to nstat-1. The equation used to
 * access the index of the associated element comes out to:
 *  bmap_index  =  tb - 1  + ta/2 (2.nstat - 3 - ta)
 */
  ib = tb - 1 + (ta * (2*nstat - 3 - ta)) / 2;
/*
 * Return the index of the sub-array descriptor.
 */
  return &an->bmap[ib];
}

/*.......................................................................
 * Read the antenna tables associated with the observation. Map each
 * table with its associated a ob->sub[] sub-array descriptor index in 
 * fob->antab[].sub. Then within each sub-array descriptor initialize
 * the associated station and [used] baseline descriptors. The
 * fob->bmap[] descriptor will also be allocated and intialized to
 * map baseline random parameters to ob->sub[].base[] baseline
 * descriptor indexes via future calls to loc_base().
 *
 * Input:
 *  ob Observation *  The descriptor of the Observation being read.
 *  fob      Fitob *  The FITS/Obsevation intermediary descriptor.
 *                    On return fob->antab[].sub will be assigned.
 * Output: 
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static int get_subarrays(Observation *ob, Fitob *fob)
{
  Antab *an;        /* Antenna table descriptor */
  int extver;       /* Antenna table extension version number */
  int isub=0;       /* Index of sub-array in ob->sub[] */
/*
 * Process each AN table.
 */
  an = fob->antab;
  for(extver=0; extver<fob->maxan; extver++,an++) {
    if(an->thdu && an->ntime>0) {
      Subarray *sub = &ob->sub[isub++];
/*
 * Record the sub-array descriptor.
 */
      fob->antab[extver].sub = sub;
/*
 * Initialize internal sub-array descriptor.
 */
      if(ini_Subarray(sub, ob->nif, an->nbase, an->nstat, an->ntime))
	return 1;
/*
 * Initialize the station descriptors in ob->sub[isub].
 */
      if(get_stations(ob, fob->fits, an))
	return 1;
/*
 * Initialize the corresponding used baseline descriptors.
 */
      if(get_baselines(ob, fob, an))
	return 1;
    };
  };
  return 0;
}

/*.......................................................................
 * Locate and index all antenna tables. Unfortunately, AIPS quite happily
 * allows both ASCII and binary versions of an antenna table with the same
 * version numbers to exist in the same file. This function compiles a list
 * of the last AN table of each version number in the file. These are
 * presumeably the most up to date versions of the tables.
 *
 * Input:
 *  fob   Fitob *  The FITS/Observation intermediary descriptor.
 *  keepant int    If true, mark all antennas in the antenna table,
 *                 plus the resulting baselines as to be kept. If
 *                 false, arrange for only those antennas and baselines
 *                 that have visibilities associated with them, to
 *                 be kept.
 * Output:
 *  fob->antab     This array of maps from AN table version number to
 *                 sub-array descriptor and AN table descriptor will be
 *                 allocated and initialized by this function.
 *                 The sub-array pointers will remain NULL until further
 *                 initialized by a call to get_subarrays() after the
 *                 observation descriptor has been allocated.
 *  fob->maxan     The dimension of fob->antab[] is maxan+1.
 *  return  int    0 - OK.
 *                 1 - Error.
 */
static int find_subarrays(Fitob *fob, int keepant)
{
  Antab *an;      /* Container of antenna table details */
  Thdu *prev;     /* The preceding antenna-table descriptor */
  Thdu *thdu;     /* The antenna table descriptor. */  
  int extver;     /* Antenna table extension version number */
  int i;
/*
 * Look at the version numbers of all AN tables to find the highest.
 */
  fob->maxan = 0;
  prev = NULL;
  while((thdu=find_table(fob->fits, "AIPS AN", -1, (Hdu *) prev)) != NULL) {
    prev = thdu;
/*
 * Ascertain the extension version number.
 */
    extver = thdu->extver;
    if(extver <= 0) {
      lprintf(stderr, "find_subarrays: Illegal AN table version number: %d\n",
	      extver);
      return 1;
    };
/*
 * Record the highest version found.
 */
    if(extver > fob->maxan)
      fob->maxan = extver;
  };
/*
 * No antenna tables found?
 */
  if(prev == NULL) {
    lprintf(stderr, "find_subarrays: No Antenna tables found.\n");
    return 1;
  };
/*
 * Allocate the map array.
 */
  fob->antab = malloc(fob->maxan * sizeof(*fob->antab));
  if(fob->antab==NULL) {
    lprintf(stderr, "find_subarrays: Insufficient memory to index AN tables\n");
    return 1;
  };
/*
 * Initialize all members of the map array to NULL.
 */
  an = fob->antab;
  for(extver=0; extver<fob->maxan; extver++,an++) {
    an->sub = NULL;
    an->thdu = NULL;
    an->integ = NULL;
    an->antrow = NULL;
    an->bmap = NULL;
    an->smap = NULL;
    an->ntime = 0;
    an->nrow  = 0;
    an->nsmax = 0;
    an->nstat = 0;
    an->nbmax = 0;
    an->nbase = 0;
    an->datutc = 0.0;
  };
/*
 * Now map AN table version numbers to AN table HDU descriptors, such that
 * the last table of a given version number is mapped.
 */
  prev = NULL;
  while((thdu=find_table(fob->fits, "AIPS AN", -1, (Hdu *) prev)) != NULL) {
    prev = thdu;
    extver = thdu->extver - 1;
    fob->antab[extver].thdu = thdu;
  };
/*
 * Ready the entries of found AN tables.
 */
  an = fob->antab;
  for(extver=0; extver<fob->maxan; extver++,an++) {
    if(an->thdu) {
/*
 * How many antennas are there?
 */
      an->nrow = numrow(an->thdu);
/*
 * Allocate and initialize the array of antenna-number -> AN table row
 * indexes in an->antrow[]. Also initialize an->nsmax.
 */
      if(get_antrow(fob->fits, an->thdu, an))
	return 1;
/*
 * Allocate and initialize the array of baseline usage maps.
 */
      an->nbmax = an->nrow * (an->nrow-1) / 2;
      an->bmap = (Anmap *) malloc(sizeof(Anmap) * an->nbmax);
      if(an->bmap==NULL) {
	lprintf(stderr,
		"find_subarrays: Insufficient memory to index baselines.\n");
	return 1;
      };
      for(i=0; i<an->nbmax; i++) {
	an->bmap[i].used = keepant;
	an->bmap[i].slot = 0;
      };
/*
 * Allocate and initialize the array of antenna usage maps.
 */
      an->smap = (Anmap *) malloc(sizeof(Anmap) * an->nrow);
      if(an->smap==NULL) {
	lprintf(stderr,
		"find_subarrays: Insufficient memory to index antennas.\n");
	return 1;
      };
      for(i=0; i<an->nrow; i++) {
	an->smap[i].used = keepant;
	an->smap[i].slot = 0;
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Read the antenna numbers of a given AIPS AN antenna table, check
 * that they are in the legal 0..255 range (limited by the encoding of
 * the BASELINE random parameter), record the max antenna number in
 * an->nmax, allocate an->antrow[] and initialize it with the row
 * indexes that correspond to each antenna.
 *
 * Input:
 *  fits       Fits * The descriptor of the FITS file.
 *  thdu       Thdu * The FITS table descriptor of the antenna table.
 *  an        Antab * The intermediary antenna table descriptor.
 * Output:
 *  return      int   0 - OK.
 *                    1 - Error.
 */
static int get_antrow(Fits *fits, Thdu *thdu, Antab *an)
{
  char *column_name;        /* The antenna-number column name */
  int icol;                 /* The antenna-number column index */
  int irow;                 /* The table row being processed */
  int i;
/*
 * Get the name used to refer to the antenna number column in the
 * type of table that the current "AIPS AN" table is using.
 */
  switch(thdu->type) {
  case F_TABLE:
    column_name = "ANT NO.";
    break;
  case F_BINTAB:
    column_name = "NOSTA";
    break;
  default:
    lprintf(stderr, "get_antrow: AN table has unknown table type.\n");
    return 1;
  };
/*
 * Search for the antenna-number column.
 */
  icol = find_column(thdu, column_name, 0);
  if(icol==0) {
    lprintf(stderr, "get_antrow: Missing %s field in AN table.\n", column_name);
    return 1;
  };
/*
 * Read the antenna number of each row and accumulate a record of the
 * max antenna number in the table.
 */
  an->nsmax = -1;
  for(irow=1; irow<=an->nrow; irow++) {
    int iant = -1;
    if(rcolumn(fits, thdu, icol, irow, DAT_INT, 1, NULL, 0, 1, &iant) != 1)
      return 1;
    if(iant < 0 || iant > 255) {
      lprintf(stderr, "get_antrow: AN-table antenna number (%d) outside legal 0..255 range.\n", iant);
      return 1;
    };
/*
 * Update the max antenna number.
 */
    if(iant > an->nsmax)
      an->nsmax = iant;
  };
/*
 * No antennas?
 */
  if(an->nsmax < 0) {
    lprintf(stderr, "get_antrow: No antennas in AN table.\n");
    return 1;
  };
/*
 * Allocate an array to index row indexes from antenna numbers.
 */
  an->antrow = (int *) malloc(sizeof(int) * (an->nsmax+1));
  if(an->antrow == NULL) {
    lprintf(stderr, "get_antrow: Insufficient memory to index AN table.\n");
    return 1;
  };
/*
 * Pre-initialize to illegal indexes so that we can detect illegal
 * antenna numbers in the data.
 */
  for(i=0; i<an->nsmax+1; i++)
    an->antrow[i] = -1;
/*
 * Read the antenna numbers again and record their 0-relative row
 * indexes in an->antrow[].
 */
  for(irow=1; irow<=an->nrow; irow++) {
    int iant;
    if(rcolumn(fits, thdu, icol, irow, DAT_INT, 1, NULL, 0, 1, &iant) != 1)
      return 1;
    an->antrow[iant] = irow-1;
  };  
  return 0;
}

/*.......................................................................
 * Search for the required and optional random group parameters and
 * store their indexes in gp (declared and defined at top of this file).
 * Also decode any UVW projection code and store it in fob->proj.
 *
 * Input:
 *  fob     Fitob *  The FITS/Observation intermediary descriptor.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int grp_parms(Fitob *fob)
{
  Phdu *phdu = (Phdu *) fob->fits->hdu;
  Gpar *gpar;  /* Descriptor of random parameter */
/*
 * Indexes not ready yet.
 */
  gp.ready = 0;
/*
 * Find the mandatory U,V,W coordinate random parameters.
 */
  if(uvw_parms(fob))
    return 1;
/*
 * Get the mandatory baseline and date random parameters.
 */
  if(loc_par("BASELINE", phdu, 1, 0, &gp.bas1, &gp.bas2) ||
     loc_par("DATE", phdu, 1, 0, &gp.dat1, &gp.dat2))
    return 1;
/*
 * Locate optional random parameters.
 */
  loc_par("FREQSEL", phdu, 0, 0, &gp.fq1, &gp.fq2);
  loc_par("INTTIM", phdu, 0, 0, &gp.dt1, &gp.dt2);
/*
 * Get the descriptor of the first DATE random parameter.
 */
  gpar = get_gpar(phdu, gp.dat1+1);
  if(gpar==NULL)
    return 1;
/*
 * Modify the date parameter offset so that it is converted to modified
 * julian date.
 */
  gpar->pzero -= 2400000.5;
/*
 * Indexes are now useable.
 */
  gp.ready = 1;
  return 0;
}

/*.......................................................................
 * Given a random group index, return a pointer to an internal static
 * struct containing the translated random parameters of that group.
 *
 * Input:
 *  fob     Fitob *  The FITS/Observation intermediary descriptor.
 *  group    long    The 0-relative index of the group, or -1 to set up
 *                   indexes.
 * Output:
 *  return Parval *  Pointer to internal static struct containing the
 *                   values read - or NULL on error.
 */
static Parval *read_pars(Fitob *fob, long group)
{
  static Parval pval;   /* Repository for decoded random-parameters */
  double *pars;         /* Pointer to buffer containing random parameters */
/*
 * Initialize random parameter indexes?
 */
  if(!gp.ready && grp_parms(fob))
    return NULL;
/*
 * Get the buffer into which the random parameters will be read.
 */
  pars = fob->pars;
/*
 * Read the requested group.
 */
  if(rgroup(fob->fits, (Phdu *)fob->fits->hdu, group, 0L, fob->npar, DAT_DBL,1,
	    NULL, pars) != fob->npar)
    return NULL;
/*
 * Get the U, V and W coords in there FITS file units (seconds of light
 * travel time over the projected baseline).
 */
  pval.uu = pars[gp.uu1] + (gp.uu2<0 ? 0.0 : pars[gp.uu2]);
  pval.vv = pars[gp.vv1] + (gp.vv2<0 ? 0.0 : pars[gp.vv2]);
  pval.ww = pars[gp.ww1] + (gp.ww2<0 ? 0.0 : pars[gp.ww2]);
/*
 * Get the time-stamp of the visibility (TAI modified JD [days]).
 * Note that the previous call to grp_parms() subtracted 2400000.5
 * from pzero such that Julian date is converted to modified Julian date.
 */
  pval.date = pars[gp.dat1] + (gp.dat2<0 ? 0.0 : pars[gp.dat2]);
/*
 * Get the FREQSEL random parameter.
 */
  pval.fqid = (gp.fq1 < 0 ? 1 :
	       (pars[gp.fq1] + (gp.fq2 < 0 ? 0 : pars[gp.fq2])));
/*
 * Get the INTTIM random parameter.
 */
  pval.inttim = (gp.dt1 < 0 ? 0.0 :
	       (pars[gp.dt1] + (gp.dt2 < 0 ? 0.0 : pars[gp.dt2])));
/*
 * Get the BASELINE random parameter and extract the sub-array and
 * antenna numbers.
 */
  {
    double basepar = pars[gp.bas1] + (gp.bas2<0 ? 0.0 : pars[gp.bas2]);
/*
 * Extract the two telescope indexes from the BASLINE random parameter.
 * Ignore the sub-array.
 * The three peices of information are encoded as:
 *  256 * antenna_1  +  antenna_2  +  0.01 * (subarray-1).
 * Where antenna_1 and antenna_2 are antenna numbers cited in the
 * antenna table.
 */
    int ibase = (int) basepar;
/*
 * The fractional part records the sub-array number.
 */
    pval.isub = (int) floor(100.0 * (basepar-ibase) + 0.5);
    if(pval.isub < 0) pval.isub = 0;
/*
 * The integral part records 256 * ta + tb.
 */
    pval.ta = (ibase >> 8);
    pval.tb = ibase - (pval.ta << 8);
  };
/*
 * Return the intialized parameters.
 */
  return &pval;
}

/*.......................................................................
 * Construct and initialize an integration grid iterator by binning
 * time-ordered groups in a UV FITS file into integration bins.
 *
 * The returned iterator is initialized such that repeated calls to
 * nxt_Intbin() will return the descriptor of each integration bin in time
 * order, and for each integration bin nxt_group() will return the index
 * of each contained group.
 *
 * Input:
 *  fob      Fitob *   The FITS/Observation intermediary descriptor.
 *  binwid  double     The bin width (seconds).
 * Output:
 *  return Intlist *   The integration list container, or NULL on error.
 */
static Intlist *bin_uvdata(Fitob *fob, double binwid)
{
  Antab *an;     /* Pointer into fob->antab[] */
  Intlist *ilist;/* The integration list container. */
  Anmap *bmap;   /* Pointer to fob->bmap[] baseline-mapping entry */
  long group;    /* Group number */
  long gcount;   /* Group count */
  long ngroup=0L;/* The minimum number of groups for complete sampling. */
  int i;
/*
 * Make sure that users are aware of what is being done to their
 * data.
 */
  if(binwid > 0.0)
    lprintf(stdout, "Binning data onto a %g second integration grid.\n",
	    binwid);
/*
 * Create the integration bin iterator.
 */
  ilist = new_Intlist(fob->maxan, fob->date.ut, binwid);
  if(ilist==NULL)
    return NULL;
/*
 * Determine the number of groups in the file.
 */
  gcount = fob->fits->hdu->gcount;
/*
 * Read the UV FITS file in group order.
 */
  for(group=fob->start_group; group<gcount; group++) {
/*
 * Read the random group parameters of the new group.
 */
    Parval *pval = read_pars(fob, group);
    if(pval==NULL)
      return del_Intlist(ilist);
/*
 * Get the baseline mapping entry for the baseline and sub-array of
 * the new group.
 */
    if((bmap=loc_base(fob, pval)) != NULL) {
/*
 * Get the offset of the times of this subarray wrt UTC.
 * Note that pval->isub has already been bounds checked by loc_base().
 */
      double datutc = fob->antab[pval->isub].datutc;
/*
 * Convert the date of the group to the number of seconds since the
 * start of the year of the observation.
 */
      double newut = (pval->date - fob->date.utc_ref) * daysec - datutc;
/*
 * Record the fact that the baseline of the new group, in the
 * associated sub-array, is sampled.
 */
      bmap->used = 1;
/*
 * Append a record of the group index to the appropriate integration bin.
 */
      if(add_group(ilist, newut, group, pval->isub))
	return del_Intlist(ilist);
    };
  };
/*
 * Count the number of baselines that are now flagged as used in
 * each sub-array. Also record the maximum of these counts,
 * count the number of sub-arrays that have associated integrations and
 * count the number of integrations, both per sub-array and as a whole.
 */
  fob->nbmax = 0;
  fob->nsub = 0;
  an = fob->antab;
  for(i=0; i<fob->maxan; i++,an++) {
/*
 * AN table i exists?
 */
    if(an->thdu) {
      int ia,ib;   /* Indexes of stations in original AN table */
      int base;    /* The index of a baseline in the original AN table */
      Anmap *bmap; /* Pointer to an->bmap */
      Anmap *smap; /* Pointer to an->smap */
/*
 * Count the number of used baselines and record their new indexes.
 */
      an->nbase = 0;
      bmap = an->bmap;
      for(base=0; base<an->nbmax; base++,bmap++) {
	if(bmap->used)
	  bmap->slot = an->nbase++;
      };
/*
 * Get the baseline and station mapping arrays.
 */
      bmap = an->bmap;
      smap = an->smap;
/*
 * Loop through the baselines and mark used stations in smap[].
 */
      for(ia=0; ia<an->nrow; ia++) {
	for(ib=ia+1; ib<an->nrow; ib++,bmap++) {
	  if(bmap->used) {
	    smap[ia].used = 1;
	    smap[ib].used = 1;
	  };
	};
      };
/*
 * Count the number of stations marked for use, and work out their
 * output indexes.
 */
      an->nstat = 0;
      smap = an->smap;
      for(ia=0; ia<an->nrow; ia++,smap++) {
	if(smap->used)
	  smap->slot = an->nstat++;
      };
/*
 * Update the max number of baselines per sub-array.
 */
      if(an->nbase > fob->nbmax)
	fob->nbmax = an->nbase;
/*
 * Determine the number of integration bins associated with the
 * antenna table.
 */
      an->ntime = ibin_count(ilist, i);
/*
 * Count used sub-arrays.
 */
      if(an->ntime > 0) {
	fob->nsub++;
/*
 * Report the results.
 */
	lprintf(stdout,
	    "AN table %d: %d integrations on %d of %d possible baselines.\n",
	     i+1, an->ntime, an->nbase, an->nbmax);
/*
 * Accumulate the overall sum of integrations in the observation.
 */
	fob->ntime += an->ntime;
/*
 * Also accumulate the number of groups that one would need to sample
 * all baselines on all integrations.
 */
	ngroup += an->ntime * an->nbase;
      } else {
	lprintf(stdout, "AN table %d: Unused.\n", i+1);
      };
    };
  };
/*
 * Abort if there are no integrations in the FITS file.
 */
  if(fob->ntime < 1) {
    lprintf(stderr, "There appear not to be any visibilities in the file.\n");
    return del_Intlist(ilist);
  };
/*
 * Produce a warning if the sampling is poor.
 */
  lprintf(stdout,
	  "Apparent sampling: %g visibilities/baseline/integration-bin.\n",
	  (double) gcount / (double) ngroup);
  if(gcount < 0.5 * ngroup)
    lprintf(stdout,
    "*** This seems a bit low - see \"help observe\" on the binwid argument.\n");
/*
 * Return the initialized iterator.
 */
  return ilist;
}

/*.......................................................................
 * Copy visibilities from a UV FITS file to the output uvdata.scr file
 * and into memory. bin_uvdata() must be called before this function,
 * to determine how to bin the input visibilities into integrations etc..
 *
 * Input:
 *  ob  Observation *  The descriptor of the Observation being initialized.
 *  fob       Fitob *  The FITS/Observation intermediary descriptor.
 * Output:
 *  return      int    Return value returned via iret argument of
 *                     uvretfn(av,iret).
 *                      0 - OK.
 *                      1 - Error.
 */
static int get_uvdata(Observation *ob, Fitob *fob)
{
  Visaver *av=NULL;    /* Visibility averager descriptor */
  Antab *an;           /* Pointer into fob->antab[] */
  Integration *integ;  /* Pointer to the latest integration */
  Intbin *ibin;        /* The integration bin being processed */
  Dpage *dp;           /* The descriptor of the output uvdata scratch file */
  double *datbuf;      /* Pointer to group data buffer */
  long irec=0;         /* Record number in scratch file */
  long igroup;         /* The index of the group being read. */
  int i;
/*
 * Get a local copy of the output uvdata file paging descriptor.
 */
  dp = ob->dp;
/*
 * Initialize the output integration record write range to encompass
 * a whole integration.
 */
  if(dp_crange(dp, 0, ob->nchan-1) || dp_irange(dp, 0, ob->nif-1) ||
     dp_brange(dp, 0, ob->nbmax-1) || dp_srange(dp, 0, ob->npol-1))
    return uvretfn(av, 1);
/*
 * Keep the user informed.
 */
  lprintf(stdout, "\nReading %ld visibilities.\n",
	  (long) fob->fits->hdu->gcount * ob->nchan * ob->nif * ob->npol);
/*
 * Get pointers into the sub-array integration arrays.
 */
  an = fob->antab;
  for(i=0; i<fob->maxan; i++,an++) {
    if(an->sub)
      an->integ = an->sub->integ;
  };
/*
 * Construct the integration bin visibility averager.
 */
  av = new_Visaver(dp, fob->binwid, fob->scatter);
  if(av==NULL)
    return uvretfn(av, 1);
/*
 * Get a pointer to the group data buffer.
 */
  datbuf = fob->data;
/*
 * Use the integration bin iterator to iterate through the groups that
 * are to be combined into integrations.
 */
  while( (ibin=nxt_Intbin(fob->ilist)) != NULL) {
/*
 * Get the descriptor of the antenna table to which the latest integration
 * belongs.
 */
    an = &fob->antab[ibin->isub];
/*
 * Get the descriptor of the integration being read.
 */
    integ = an->integ++;
/*
 * Initialize the averaged output record.
 */
    if(av_newint(av, integ->vis, integ->sub->nbase, irec))
      return uvretfn(av, 1);
/*
 * Record the integration time-stamp. This is the number of seconds into
 * the year in which the observation started.
 */
    integ->ut = ibin->ut;
/*
 * Also record the output uvdata record number.
 */
    integ->irec = irec++;
/*
 * Read each group to be binned into the output integration.
 */
    while( (igroup=nxt_group(ibin)) != -1) {
      Anmap *bmap;        /* Baseline index mapping container */
/*
 * Read the random-group parameters of the new group.
 */
      Parval *pval = read_pars(fob, igroup);
      if(pval==NULL)
	return uvretfn(av, 1);
/*
 * Does the new group contain data for the sub-array of the current
 * integration bin, and cite a useable baseline?
 */
      if(pval->isub == ibin->isub && (bmap=loc_base(fob, pval))!=NULL) {
	int base;           /* The ob index of the baseline being read */
	int xif;            /* Index of IF */
	int ch;             /* Index of spectral-line channel */
	int pol;            /* Index of polarization */
	float group_wt=0.0f;/* Overall weight of the current group */
/*
 * Use it to get the baseline index in the observation structure
 * and output scratch file, to be read into.
 */
	base = bmap->slot;
/*
 * Read the visibility data of the new group.
 */
	if(get_data(fob, igroup)==NULL)
	  return uvretfn(av, 1);
/*
 * Add to the weighted running mean of the visibilities in the output
 * buffer for baseline 'base' at each spectral-line channel and polarization
 * in each IF.
 */
	for(xif=0; xif<ob->nif; xif++) {  /* Loop through the IFs */
	  Dif *ifs = &dp->ifs[xif];
	  long ifpos = xif * ax.iinc;
	  for(ch=0; ch<ob->nchan; ch++) { /* Spectral-line channels */
	    Cvis *cvis = ifs->chan[ch].base[base].pol;
	    long chpos = ifpos + ch * ax.finc;
	    for(pol=0; pol<ob->npol; pol++,cvis++) { /* Polarizations */
	      long datpos = chpos + pol * ax.sinc;
/*
 * Extract the weighted complex visibility. Apply the sign of AIPS WTSCAL
 * to the weight to instate the normal sign convention of -ve values for
 * flagged data.
 */
	      float re = datbuf[datpos];                 /* Real part */
	      float im = datbuf[datpos+1];               /* Imaginary part */
	      float wt = datbuf[datpos+2] * fob->wtsign; /* Weight */
/*
 * Ignore deleted data.
 */
	      if(wt!=0.0f) {
/*
 * Accumulate the weighted running mean visibility for the latest
 * integration bin. Initialize means on the first visibility and on
 * the first un-flagged visibility.
 */
		if(cvis->wt==0.0f || (cvis->wt<0.0f && wt>0.0f)) {
		  cvis->wt = wt;
		  cvis->im = im;
		  cvis->re = re;
/*
 * Optionally accumulate running means for use in estimates of the data scatter.
 */
		  if(av->scatsum) {
		    Scatsum *scatsum = &av->scatsum[cvis-dp->cvis];
		    scatsum->sqr_mean = (re*re+im*im);
		    scatsum->nsum = 1;
		  };
		}
/*
 * Accumulate the running mean U,V,W, using flagged visibilities only until
 * the flagged means are reset (see above) to un-flagged by the first good
 * visibility.
 */
		else if(wt > 0.0f || cvis->wt < 0.0f) {
		  float runwt = wt / (cvis->wt += wt);
		  cvis->im += runwt * (im - cvis->im);
		  cvis->re += runwt * (re - cvis->re);
/*
 * Optionally accumulate running means for use in estimates of the data scatter.
 */
		  if(av->scatsum) {
		    Scatsum *scatsum = &av->scatsum[cvis-dp->cvis];
		    scatsum->sqr_mean += runwt * (re*re+im*im - scatsum->sqr_mean);
		    scatsum->nsum++;
		  };
		};
/*
 * Accumulate the overall weight of the group.
 */
		if(group_wt==0.0f || (group_wt<0.0f && wt>0.0f))
		  group_wt = wt;
		else if(wt>0.0f || group_wt<0.0f)
		  group_wt += wt;
	      };
	    };
	  };
	};
/*
 * Accumulate the weighted mean U,V,W coordinates of the current baseline.
 */
	av_uvwt(av, pval->uu, pval->vv, pval->ww, group_wt, pval->inttim, base);
      };
    };
/*
 * Set U,V and W to zero on baselines which were totally un-sampled, and
 * if requested, replace output weights with those deduced from the data
 * scatter.
 */
    if(av_endint(av))
      return uvretfn(av, 1);
/*
 * Write the accumulated output integration record to ob->dp.
 */
    if(dp_write(dp, integ->irec))
      return uvretfn(av, 1);
  };
  return uvretfn(av, 0);
}

/*.......................................................................
 * Private cleanup function of get_uvdata().
 */
static int uvretfn(Visaver *av, int iret)
{
  av = del_Visaver(av);
  return iret;
}

/*.......................................................................
 * Locate the U,V and W coordinate random parameters, record their
 * locations in gp.* and decode and record their projection type in
 * fob->proj. Also if UU-L is found instead of UU, fix it by scaling
 * poff and pscal by the CRVAL value of the FREQ axis.
 *
 * Input:
 *  fob     Fitob *  The FITS/Observation intermediary descriptor.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int uvw_parms(Fitob *fob)
{
  char uvwname[9];     /* Full space padded 8 character parameter name */
  Phdu *phdu;          /* The Primary HDU descriptor */
  double xfreq;        /* The frequency needed to re-cale UU-L */
  char *cptr;          /* Pointer into a random parameter name */
  int i;
/*
 * List the attributes of each parameter to be acquired, in the order, U,V,W.
 */
  typedef struct {
    char *name;      /* The main part of the name (UU,VV,WW) */
    int nfound;      /* The number of matching parameters */
    int ip[2];       /* The 1-relative indexes of the parameters */
    double scale[2]; /* Any scale factor required to convert to light-seconds */
  } UVWpar;
  UVWpar uvwpar[3] = {
    {"UU", 0, {0,0}, {0.0,0.0}},
    {"VV", 0, {0,0}, {0.0,0.0}},
    {"WW", 0, {0,0}, {0.0,0.0}}
  };
/*
 * Record the UVW projection.
 */
  Proj uvwproj = PRJ_NON;
/*
 * Get the primary HDU descriptor.
 */
  phdu = (Phdu *) fob->fits->hdu;
/*
 * Get the reference value of the frequency axis in case we have to
 * scale UU-L etc.. to UU etc..
 */
  {
    Imaxis *axis;
    int fpos = find_axis(phdu, "FREQ", 0, 1);
    if(fpos <= 0) {
      lprintf(stderr, "uvw_parms: Unable to find the FREQ axis.\n");
      return 1;
    };
    if((axis=get_axis(phdu, fpos))==NULL)
      return 1;
    xfreq = axis->crval;
  };
/*
 * Attempt to locate each of UU,VV and WW.
 */
  for(i=0; i<3; i++) {
    UVWpar *par = uvwpar + i;
    int ip;
/*
 * Check each random parameter name against par->name.
 */
    for(ip=0; ip<phdu->pcount; ip++) {
      char *name=gpar_name(phdu, ip+1);
/*
 * par->name is actually just a 2 letter prefix of a number of possible
 * matching names.
 */
      if(name && strncmp(name, par->name, 2)==0) {
	int found = 0;      /* True if the parameter name is acceptable */
	double scale = 1.0; /* The UVW coordinate scale factor */
/*
 * Convert the name into a right justified, 8 character, space padded string.
 * UVW coordinate parameter names are encoded as a 4 character type, optionally
 * followed by a 4 character projection. Padding to the maximum of 8
 * characters thus allows us to ignore the possibility of smaller lengths.
 */
	sprintf(uvwname, "%-8.8s", name);
/*
 * Convert hyphens to spaces. This cuts down on the number of options that
 * we need to check.
 */
	for(cptr=uvwname; *cptr; cptr++) {
	  if(*cptr=='-')
	    *cptr = ' ';
	};
/*
 * Check for parameter names of the form "UU  *".
 */
	if(strncmp(uvwname+2, "  ",  2)==0) {
	  found = 1;
/*
 * Check for parameter names of the form "UU-L*" (ie. "UU L*").
 */
	} else if(strncmp(uvwname+2, " L", 2)==0) {
	  found = 1;
	  scale = 1.0/xfreq;  /* Convert from wavelengths to light-seconds */
	};
/*
 * Did we get a new match?
 */
	if(found) {
/*
 * Too many matches?
 */
	  if(++par->nfound > 2) {
	    lprintf(stderr, "Too many (>2) %s random parameters.\n", par->name);
	    return 1;
	  };
/*
 * Record the new parameter.
 */
	  par->ip[par->nfound-1] = ip+1;
	  par->scale[par->nfound-1] = scale;
/*
 * Skip white-space past character 4, up to the first character of the
 * projection type name.
 */
	  for(cptr = uvwname + 4; *cptr && isspace((int)*cptr); cptr++)
	    ;
	  {
	    Proj proj = *cptr ? name_Proj(cptr) : PRJ_SIN;
	    if(proj == PRJ_NON) {
	      lprintf(stderr,
	      "uvw_parms: %c coordinate projection \"%s\" is not recognized.\n",
		      par->name[0], cptr);
	      return 1;
	    };
/*
 * Record the first projection found and check subsequent projection types
 * against it.
 */
	    if(uvwproj==PRJ_NON) {
	      uvwproj = proj;
	    } else if(uvwproj != proj) {
	      lprintf(stderr,
	     "uvw_parms: Inconsistent UVW coordinate projections, %s and %s.\n",
		      Proj_name(uvwproj), Proj_name(proj));
	      return 1;
	    };
	  };
	};
      };
    };
/*
 * Parameter not found?
 */
    if(!par->nfound) {
      lprintf(stderr, "Failed to find %s random parameter.\n", par->name);
      return 1;
    };
/*
 * Convert units to light-seconds.
 */
    {
      int pp;
      for(pp=0; pp<par->nfound; pp++) {
	Gpar *gpar = get_gpar(phdu, par->ip[pp]);
	if(!gpar)
	  return 1;
	gpar->pscal *= par->scale[pp];
	gpar->pzero *= par->scale[pp];
      };
    };
  };
/*
 * Store the parameter indexes for later use.
 */
  {
    UVWpar *par = uvwpar;
    gp.uu1 = par->ip[0] - 1;
    gp.uu2 = par->ip[1] - 1;
    par++;
    gp.vv1 = par->ip[0] - 1;
    gp.vv2 = par->ip[1] - 1;
    par++;
    gp.ww1 = par->ip[0] - 1;
    gp.ww2 = par->ip[1] - 1;
/*
 * Record the coordinate projection.
 */
    fob->proj = uvwproj;
  };
  return 0;
}

/*.......................................................................
 * Read the reference antenna and R-L phase difference header keywords
 * from an antenna table.
 *
 * Input:
 *  fits       Fits *  The descriptor of the FITS file.
 *  hdu         Hdu *  The Hdu descriptor of the antenna table.
 *  sub    Subarray *  The sub-array associated with the antenna table.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
static int rd_p_refant(Fits *fits, Hdu *hdu, Subarray *sub)
{
  Fitkey key;        /* Descriptor of FITS keyword/value pair */
/*
 * See if there is a P_REFANT header keyword. This specifies a reference
 * antenna used in polarization calibration. Associated with this are
 * sub->nif keywords named P_DIFFnn where nn=1..nif.
 */
  if(get_key(fits, hdu, "P_REFANT", DAT_INT, LOOP_SEEK,
	     &key)==KEY_FOUND) {
    enum {AN_P_DIFF=0};
    static Fitkey p_keys[] = {
      {"P_DIFF", AN_P_DIFF, 0, DAT_DBL, NULL, NULL}
    };
/*
 * Record the reference antenna number.
 */
    sub->p_refant = KEYINT(key);
/*
 * Search for the associated P_DIFFnn keywords.
 */
    new_hline(hdu, 0);   /* Rewind header */
    while(next_key(fits, hdu, p_keys, sizeof(p_keys)/sizeof(Fitkey),
		   EOH_SEEK, &key) == KEY_FOUND) {
      if(key.keyid==AN_P_DIFF) {
	int cif = key.extn - 1;
	if(cif >= 0 && cif < sub->nif)
	  sub->p_diff[cif] = KEYDBL(key);
      };
    };
  };
  return 0;
}

/*.......................................................................
 * Return non-zero if the specified string is NULL, contains only spaces
 * or has zero length.
 *
 * Input:
 *  string    char *  The string to be checked.
 * Output:
 *  return     int    0 - The string is not empty.
 *                    1 - The string is empty.
 */
static int string_is_empty(char *string)
{
/*
 * If the string is NULL, treat it as being empty.
 */
  if(!string)
    return 1;
/*
 * Find the first non-space character in the string.
 */
  while(*string && isspace((int)*string))
    string++;
/*
 * If we reached the end of the string without finding any non-space
 * characters, then the string is empty.
 */
  return *string ? 0 : 1;
}

/*.......................................................................
 * The visibilities of each sub-array can have different time-systems.
 * We need this information before we can index the data and before
 * we can figure out the start time of the observation. This function
 * reads the time system and offset from each antenna table, or records
 * the offset in fob->antab[*].datutc.
 *
 * Input:
 *  fob      Fitob *  The FITS/Observation intermediary descriptor.
 *  iatutc  double    The default value to record for the difference
 *                    between IAT and UTC, where not provided.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int get_subarray_time_systems(Fitob *fob, double iatutc)
{
  Fitkey key;        /* The descriptor of a FITS keyword/value pair */
  int extver;        /* Antenna table extension version number */
/*
 * Create a table of the binary antenna-table keywords that are relevant
 * to date interpretation.
 */
  enum {DATUTC, IATUTC, TIMSYS};
  static Fitkey ankeys[]={
    {"DATUTC", 0, DATUTC,   DAT_DBL, NULL, NULL},
    {"IATUTC", 0, IATUTC,   DAT_DBL, NULL, NULL},
    {"TIMSYS", 0, TIMSYS,   DAT_STR, NULL, NULL},
  };
/*
 * Process each of the antenna tables that were selected in find_subarrays().
 */
  for(extver=0; extver<fob->maxan; extver++) {
    Antab *an = fob->antab + extver;
    Thdu *thdu = an->thdu;
/*
 * Is this one of the tables that was selected in find_subarrays()?
 */
    if(thdu) {
      char opt_timsys[4] = {0};  /* A TIMSYS keyword value */
      double opt_iatutc = 0;     /* A IATUTC keyword value */
      an->datutc = 0.0;          /* A DATUTC keyword value */
/*
 * ASCII tables imply IAT.
 */
      switch(thdu->type) {
      case F_TABLE:
	an->datutc = iatutc;
	break;
/*
 * In binary tables, the optional TIMSYS keyword can be used to override
 * the IAT default.
 */
      case F_BINTAB:
	new_hline((Hdu *) thdu, 0);   /* Rewind header */
/*
 * Read the keys that are listed above in the definition of ankeys[].
 */
	while(next_key(fob->fits, (Hdu *) thdu, ankeys,
		       sizeof(ankeys)/sizeof(Fitkey), EOH_SEEK, &key) == 0) {
	  switch(key.keyid) {
	  case TIMSYS:
	    stripcpy(opt_timsys, sizeof(opt_timsys),
		     KEYSTR(key), strlen(KEYSTR(key)));
	    break;
	  case IATUTC:
	    opt_iatutc = KEYDBL(key);
	    break;
	  case DATUTC:
	    an->datutc = KEYDBL(key);
	    break;
	  default:
	    break;
	  };
	};
/*
 * The default for TIMSYS is "IAT" according to Going AIPS.
 */
	if(string_is_empty(opt_timsys))
	  strcpy(opt_timsys, "IAT");
/*
 * Old files have IATUTC, and don't have DATUTC.
 * In such cases, when the time system of the data is IAT, assign
 * iatutc to the recorded datutc.
 */
	if(strcmp(opt_timsys, "IAT") == 0) {
	  if(an->datutc <= 0.0)
	    an->datutc = (opt_iatutc != 0.0) ? opt_iatutc : iatutc;
/*
 * For UTC time-systems datutc should obviously be zero.
 */
	} else if(strcmp(opt_timsys, "UTC") == 0) {
	  if(an->datutc != 0.0) {
	    lprintf(stderr,
	      "Warning: Resetting DATUTC from %g to 0, because TIMSYS='UTC'.\n",
		    an->datutc);
	    };
	    an->datutc = 0.0;
/*
 * If an unknown time system is encountered and DATUTC doesn't provide
 * its offset from UTC, complain.
 */
	} else if(an->datutc == 0) {
	  lprintf(stderr, "Warning: TIMSYS '%s' is unknown and DATUTC=0.\n",
		  opt_timsys);
	};
	break;
      default:
	lprintf(stderr,
	      "get_subarray_time_systems: AN table has unknown table type.\n");
	return 1;
      };
    };
  };
  return 0;
}

