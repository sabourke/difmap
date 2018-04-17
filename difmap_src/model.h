/*
 * Define an object type for storing vlbi model components.
 */

#ifndef model_h
#define model_h

/*
 * Define enumerators for a bitmap of free model-component parameters.
 * Each enumeration constant must be a different power of two up to a
 * maximum value of 2^15.
 */
typedef enum {
  M_FLUX=1, M_CENT=2, M_MAJOR=4, M_RATIO=8, M_PHI=16, M_SPCIND=32
} Modpar;

/*
 * Define an enumeration of recognised model-component types.
 */
typedef enum {M_DELT, M_GAUS, M_DISK, M_ELLI, M_RING, M_RECT, M_SZ} Modtyp;

typedef struct mod_type {
  Modtyp type;    /* Model component type */
  int freepar;    /* Bitmap of Modpar values used to mark free parameters */
  float flux;     /* Flux of component */
  float x;        /* RA relative coordinate of component centroid (radians) */
  float y;        /* Dec relative coordinate of component centroid (radians) */
  float major;    /* Major axis of component (radians) */
  float ratio;    /* Axial ratio (minor/major), < 1.0 */
  float phi;      /* Position angle of major axis (radians N->E) */
  float freq0;    /* The reference frequency to use with spectral indeces */
  float spcind;   /* The spectral index of the component flux */
  struct mod_type *next;
} Modcmp;

/*
 * Define a container for a list of components.
 */

typedef struct Model {
  int ncmp;        /* Number of model components in list */
  int issqd;       /* True if the contained model is currently squashed */
  int isdelt;      /* True if the model is only comprised of deltas */
  float flux;      /* Total flux in model */
  Modcmp *head;    /* Head of cmp list */
  Modcmp *tail;    /* Tail of cmp list */
} Model;


/* Model-object constructor */

Model *new_Model(void);

/* Model-object destructor */

Model *del_Model(Model *mod);

/* Clear the contents of a model */

Model *clr_Model(Model *mod);

/* Make a copy of a model */

Model *cpy_Model(Model *mod);

/* Model method functions */

Modcmp * new_cmp(int modnum);
Modcmp *del_cmp(Modcmp *cmp);

Modcmp *add_xycmp(Model *mod, int docomp, int freepar, float flux, float x,
		  float y, float major, float ratio, float phi, Modtyp type,
		  float freq0, float spcind);

Modcmp *add_cmp(Modcmp *cmp, Model *mod, int docomp);

Model *add_mod(Model *mod, Model *old, int docomp, int append);

Modcmp *rem_cmp(Model *mod, Modcmp *prev, Modcmp *cmp);

int rmodel(Model *mod, float east, float north, int docomp, char *modfile);

int wmodel(Model *mod, float east, float north, int docut, float cut, FILE *fd);

Model *cut_mod(Model *mod, float cut);

Model *squash(Model *mod);

void shiftmod(Model *mod, float east, float north);

/* Plot the fixed and/or variable components of a model */

int modplot(Model *mod, int dofix, int dovar, float xa, float xb, float ya, float yb);

/* Plot a model component if it lies in the given area */

int cmpplot(Modcmp *cmp, float xa, float xb, float ya, float yb, int erase);

/* Place the variable components of amod+bmod in bmod, and the rest in amod */

int var_mod(Model *amod, Model *bmod);

/* Allow the user to edit a model in an external editor */

Model *ed_model(Model *mod);

/* Read a single model component into a specified model */

typedef enum { /* Enumerate the return status of read_Modcmp() */
  CMP_READ,    /* A new component was read successfully. */
  CMP_EMPTY,   /* An empty line, EOF, or a component with zero flux was read. */
  CMP_ERROR    /* A fatal error occurred. */
} RModcmp;

RModcmp read_Modcmp(Model *mod, float east, float north, int docomp,
		    char *modfile, FILE *fp, int *nline);

#endif
