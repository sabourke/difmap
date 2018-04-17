#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "logio.h"
#include "vlbconst.h"
#include "model.h"
#include "scrfil.h"
#include "freelist.h"


/*
 * A freelist of model components will be allocated and assigned to
 * modcmp_memory. Whenever this runs out of memory for new models,
 * it allocates a block of MODCMP_INC more models.
 */
#define MODCMP_INC 512

static FreeList *modcmp_memory = NULL;  /* A (Modcmp *) freelist */

/*
 * A freelist of model containers will be allocated and assigned to
 * model_memory. Whenever this runs out of memory for new models,
 * it allocates a block of MODEL_INC more models.
 */
#define MODEL_INC 50

static FreeList *model_memory = NULL;   /* A (Model *) freelist */

/*.......................................................................
 * Get a new model component structure from the free-list. Also
 * initialise it so that all members are zero, including the 'next' link.
 *
 * Input:
 *   modnum  int       The number of this model component - used only
 *                     for error reporting on failure to allocate memory.
 * Output:
 *   return  Modcmp *  The new model component, or NULL on error.
 */
Modcmp *new_cmp(int modnum)
{
  Modcmp *newcmp;   /* The new Model component */
/*
 * Allocate the freelist that is used for allocating model components?
 */
  if(!modcmp_memory) {
    modcmp_memory = new_FreeList("new_cmp", sizeof(Modcmp), MODCMP_INC);
    if(!modcmp_memory)
      return NULL;
  };
/*
 * Allocate a new component from the freelist.
 */
  newcmp = new_FreeListNode("new_cmp", modcmp_memory);
  if(!newcmp)
    return NULL;
/*
 * Initialise the model component before returning it.
 */
  newcmp->next=NULL;
  newcmp->flux=0.0;
  newcmp->x=0.0;
  newcmp->y=0.0;
  newcmp->major=0.0;
  newcmp->ratio = 0.0;
  newcmp->phi = 0.0;
  newcmp->freq0 = 0.0;
  newcmp->spcind = 0.0;
  return newcmp;
}

/*.......................................................................
 * Place a model component back on the free-list. NOTE. The caller must
 * relink around this component (use rem_cmp()) before sending it to
 * this routine - otherwise the caller's model component list will
 * thereafter dive into the free-list via the modified ->next member of
 * this component.
 * eg. cmp = del_cmp(rem_cmp(mod,prev,cmp));  Where prev->next=cmp .
 *
 * Input:
 *   cmp  Modcmp *  The unlinked model component to be returned to the
 *                  free-list.
 * Output:
 *   Modcmp * 0     Allways.
 */
Modcmp *del_cmp(Modcmp *cmp)
{
  return del_FreeListNode("del_cmp", modcmp_memory, cmp);
}

/*.......................................................................
 * Dynamically allocate and initialize a vlbi Model container class
 * and return its pointer on success, or NULL on failure.
 *
 * Output:
 *   return   Model *   The initialised Model container or NULL on error.
 */
Model *new_Model(void)
{
  Model *mod;
/*
 * Create a freelist from which to allocate Model containers?
 */
  if(!model_memory) {
    model_memory = new_FreeList("new_Model", sizeof(Model), MODEL_INC);
    if(!model_memory)
      return NULL;
  };
/*
 * Allocate the model object.
 */
  mod = new_FreeListNode("new_Model", model_memory);
  if(!mod)
    return NULL;
/*
 * Initially the container is empty.
 */
  mod->ncmp = 0;
  mod->head  = 0;
  mod->tail = 0;
/*
 * Until proved otherwise assume that the model is squashed and formed
 * from delta components only.
 */
  mod->isdelt = 1;
  mod->issqd = 1;
  mod->flux = 0.0f;
/*
 * return the initialized model.
 */
  return mod;
}

/*.......................................................................
 * Delete a Model container and its contents. This involves free()ing
 * all memory allocated and returning components to the free-list.
 *
 * Input:
 *   mod  Model *   The Model to be deleted.
 * Output:
 *   return Model * The deleted model (ie. a NULL model pointer!).
 */
Model *del_Model(Model *mod)
{
/*
 * Do nothing if the Model has already been deleted.
 */
  if(mod == NULL)
    return mod;
/*
 * Return all model components to the free-list.
 */
  clr_Model(mod);
/*
 * Return the model container to the freelist.
 */
  return del_FreeListNode("del_Model", model_memory, mod);
}

/*.......................................................................
 * Delete all the components of a model, leaving an model container that
 * contains 0 components. This is an alternative to del_Model() which
 * also deletes the container.
 *
 * Input:
 *  mod    Model *  The model to be cleared (NULL is OK).
 * Output:
 *  return Model *  The cleared model container (mod).
 */
Model *clr_Model(Model *mod)
{
  Modcmp *cmp;   /* The current model component being free'd. */
  Modcmp *next;  /* The component following 'cmp' */
/*
 * Do nothing if the Model is deleted.
 */
  if(mod == NULL)
    return mod;
/*
 * Return all model components to the free-list.
 */
  for(cmp=mod->head; cmp != NULL; cmp=next) {
    next = cmp->next;
    del_cmp(cmp);
  };
/*
 * Reset the model parameters to represent its new empty state.
 */
  mod->ncmp = 0;
  mod->head  = 0;
  mod->tail = 0;
/*
 * Until proved otherwise assume that the model is squashed and formed
 * from delta components only.
 */
  mod->isdelt = 1;
  mod->issqd = 1;
  mod->flux = 0.0f;
/*
 * Return the empty container.
 */
  return mod;
}

/*.......................................................................
 * Append a new model component to a given model. (See also add_rtcmp()
 * which adds models in terms of polar coordinte positions rather
 * than x,y).
 *
 * Input:
 *   mod     Model *  The model to which the component is to be added.
 *   docomp    int    If true then merge with an existing component
 *                    if both are delta functions at the same location.
 *   freepar   int    A bitmap union of Modfree values, denoting which if
 *                    any of the model component parameters are to be
 *                    free to vary in modelfit.
 *   flux    float    The component flux at freq0.
 *   x       float    The X-position of the component (radians).
 *   y       float    The Y-position of the component (radians).
 *   major   float    Major axis of component (radians).
 *   ratio   float    Axial ratio (minor/major ie < 1).
 *   phi     float    Position angle of major axis (radians, N->E).
 *   type   Modtyp    The model type (M_DELT, M_GAUS, M_DISK, M_ELLI,
 *                    M_RING, M_RECT, M_SZ.
 *   freq0   float    The frequency corresponding to the above flux.
 *   spcind  float    The spectral index of the component flux.
 * Output:
 *   return Modcmp *  The model component added or NULL on error.
 */
Modcmp *add_xycmp(Model *mod, int docomp, int freepar, float flux, float x,
		  float y, float major, float ratio, float phi, Modtyp type,
		  float freq0, float spcind)
{
  Modcmp *cmp;     /* The new component */
/*
 * Check the reference frequency.
 */
  if((spcind != 0.0 || (freepar & M_SPCIND)) && freq0 <= 0.0) {
    lprintf(stderr,
	    "add_xycmp: Invalid model component reference frequency.\n");
    return NULL;
  };
/*
 * If the component is not a continuum delta function then the model
 * is no longer comprised only of simple deltas - record this.
 * Similarly - only if the docomp option has been selected will a
 * squashed model remain squashed.
 */
  mod->isdelt = mod->isdelt && type==M_DELT && spcind==0.0 && 
    !(freepar & M_SPCIND);
  mod->issqd  = mod->issqd && docomp;
/*
 * If the squash option is selected then first check if the new component
 * can be combined with a preceding one. If so, no new component is
 * required and simply add the new flux to the matching component.
 */
  if(docomp && type == M_DELT) {
    for(cmp=mod->head; cmp != NULL; cmp=cmp->next) {
      if(x == cmp->x && y == cmp->y && cmp->type == M_DELT &&
	 spcind == cmp->spcind) {
/*
 * If the equal spectral index of the two components is non-zero,
 * get the flux to be added at the reference frequency of the target
 * component.
 */
	if(spcind != 0.0 && freq0 != cmp->freq0)
	  flux *= pow(cmp->freq0 / freq0, spcind);
	cmp->freepar |= freepar;
	cmp->flux += flux;
	mod->flux += flux;
	return cmp;
      };
    };
  };
/*
 * Get a new component from the free list.
 */
  cmp = new_cmp(mod->ncmp+1);
  if(cmp == NULL)
    return cmp;
/*
 * Fill in the component attributes.
 */
  cmp->type=type;
  cmp->freepar = freepar;
  cmp->x = x;
  cmp->y = y;
  cmp->flux = flux;
  cmp->major=major;
  cmp->ratio=ratio;
  cmp->phi=phi;
  cmp->spcind = spcind;
  cmp->freq0 = freq0;
/*
 * Insert the new component and return it.
 */
  return add_cmp(cmp,mod,0);
}

/*.......................................................................
 * Read the model components file 'modfile' and insert its components
 * into Model 'mod'.
 *
 * Input:
 *  mod    Model *  An exisiting Model into which the new components
 *                  should be inserted.
 *  east   float    Any eastward X-axis offset to be added while reading.
 *  north  float    Any northward Y-axis offset to be added while reading.
 *  docomp   int    If true then merge delta components that have the
 *                  same location.
 *  modfile char *  The name of the model file.
 * Output:
 *   return   int   0 - OK.
 *                  1 - Error.
 */
int rmodel(Model *mod, float east, float north, int docomp, char *modfile)
{
  FILE *fd;       /* Model-file pointer */
  int nread=0;    /* Number of components successfully read from file */
  int nline;      /* Current line number in modfile. */
/*
 * Open the model file.
 */
  fd = fopen(modfile, "r");
  if(fd == NULL) {
    lprintf(stderr, "Unable to open model file: %s\n", modfile);
    return 1;
  };
/*
 * Each component occupies a single line in the file.
 */
  for(nline=1; !feof(fd); nline++) {
    switch(read_Modcmp(mod, east, north, docomp, modfile, fd, &nline)) {
    case CMP_READ:
      nread++;
      break;
    case CMP_EMPTY:
      break;
    case CMP_ERROR:
      fclose(fd);
      return 1;
      break;
    };
  };
/*
 * Succesfull finish report the number of components read.
 */
  lprintf(stdout, "A total of %d model components were read from file %s\n",
	  nread, modfile);
/*
 * Close model file.
 */
  fclose(fd);
  return 0;
}

/*.......................................................................
 * Read a single line of a model file, and add the resulting component
 * to the specified model.
 *
 * Input:
 *  mod      Model *  The model being read.
 *  east     float    Any eastward X-axis offset to be added while reading.
 *  north    float    Any northward Y-axis offset to be added while reading.
 *  docomp     int    If true then merge delta components that have the
 *                    same location.
 *  modfile   char *  The name of the input file.
 *  fp        FILE *  The input stream connected to the input file.
 * Input/Output:
 *  nline      int *  The line number will be maintained in *nline. If the
 *                    caller does any reading it is expected to increment
 *                    this whenever a newline is seen.
 *  
 * Output:
 *  return RModcmp    One of the values:
 *                     CMP_READ  -  A new component was read successfully.
 *                     CMP_EMPTY -  An empty line, EOF, or a component with
 *                                  zero flux was read.
 *                     CMP_ERROR -  A fatal error occurred.
 */
RModcmp read_Modcmp(Model *mod, float east, float north, int docomp,
		    char *modfile, FILE *fp, int *nline)
{
  int nfield=0;   /* The number of fields read from the line */
  int was_eol=0;  /* Set to true at end of the line */
  int c;          /* Character returned by fgetc() */
/*
 * Assign defaults to the attributes of the next component.
 */
  float flux = 0.0f;    /* Flux density (map units) */
  float radius = 0.0f;  /* Radius (mas) */
  float theta = 0.0f;   /* Position angle (degrees) North through East */
  float major = 0.0f;   /* Major axis (mas) */
  float ratio = 1.0f;   /* Axial ratio (minor/major) <= 1.0 */
  float phi = 0.0f;     /* Angle of major axis (degrees) North though East */
  float freq0 = 0.0f;   /* The reference frequency of the component */
  float spcind = 0.0f;  /* The spectral index of the component flux */
  int type = 0;         /* Component type */
  int freepar=0;        /* Bitmap union of free-parameter type enumeration */
/*
 * Read up to 9 fields from the next line.
 */
  while(!was_eol) {
    int m_type = 0;  /* Free parameter type enumeration */
/*
 * Skip intermediate spaces and escaped newlines.
 */
    do {
      c = fgetc(fp);
/*
 * Did we encounter an escape character?
 */
      while(c == '\\') {
/*
 * Search ahead for the next non-space character.
 */
	do {
	  c = fgetc(fp);
	} while(c==' ' || c=='\t');
/*
 * If it is a newline, skip it and increment the line counter.
 */
	if(c=='\n') {
	  c = fgetc(fp);
	  (*nline)++;
	} else {
	  lprintf(stderr, "rmodel: Unespected escape character on line %d\n",
		  *nline);
	  return CMP_ERROR;
	};
      };
    } while(c==' ' || c=='\t');
/*
 * Abort further reads if the end of line was reached.
 */
    if(c==EOF || c=='\n') {
      was_eol = 1;
    } else {
      ungetc(c, fp);
/*
 * Read the next field.
 */
      switch(nfield) {
      case 0:
	was_eol = fscanf(fp, "%f", &flux) != 1;
	m_type = M_FLUX;
	break;
      case 1:
	was_eol = fscanf(fp, "%f", &radius) != 1;
	m_type = M_CENT;
	break;
      case 2:
	was_eol = fscanf(fp, "%f", &theta) != 1;
	m_type = M_CENT;
	break;
      case 3:
	was_eol = fscanf(fp, "%f", &major) != 1;
	m_type = M_MAJOR;
	break;
      case 4:
	was_eol = fscanf(fp, "%f", &ratio) != 1;
	m_type = M_RATIO;
	break;
      case 5:
	was_eol = fscanf(fp, "%f", &phi) != 1;
	m_type = M_PHI;
	break;
      case 6:
	was_eol = fscanf(fp, "%d", &type) != 1;
	break;
      case 7:
	was_eol = fscanf(fp, "%f", &freq0) != 1;
	break;
      case 8:
	was_eol = fscanf(fp, "%f", &spcind) != 1;
	m_type = M_SPCIND;
	break;
      default:
	lprintf(stderr, "rmodel: Too many fields on line %d\n", *nline);
	return CMP_ERROR;
	break;
      };
/*
 * End of line?
 */
      if(was_eol) {
/*
 * Was the line ended by a comment, or by an error?
 */
	if(c != '!') {
	  lprintf(stderr,"rmodel: Error at field %d on line %d of file: %s\n",
		  nfield+1, *nline, modfile);
	  return CMP_ERROR;
	};
/*
 * Skip to the start of the next line.
 */
	while((c=fgetc(fp)) != '\n' && c!=EOF);
/*
 * A new field has been succesfully read.
 */
      } else {
/*
 * Accumulate the count of the number of fields sucessfully read.
 */
	nfield++;
/*
 * Fields that are post-fixed with a 'v' or 'V', are to be marked as
 * free-parameters for model fitting. Accumulate the bitmap union of
 * free parameter types.
 */
	c=getc(fp);
	if(c=='v' || c=='V')
	  freepar |= m_type;
	else if(c!=EOF)
	  ungetc(c,fp);
      };
    };
  };
/*
 * Act upon the different possible number of fields read.
 */
  if(nfield>0 && flux != 0.0f) {
/*
 * If the model type was ommitted, determine it from the rules given in
 * the "Introduction to Caltech VLBI programs" document.
 */
    if(nfield<7) {
      if(nfield <= 3 || major == 0.0)
	type = 0;		/* DELTA function component */
      else
	type = 1;		/* GAUSSIAN component */
    };
/*
 * If the major axis is zero then convert the component to a delta
 * component.
 */
    if(major==0.0)
      type = 0;
/*
 * Discard irrelevant parameters from delta components.
 */
    if(type==0) {
      major = 0.0f;
      ratio = 1.0f;
      phi = 0.0f;
      freepar &= ~(M_MAJOR | M_RATIO | M_PHI);
    };
/*
 * Check the legality of the component type.
 */
    if(type < M_DELT || type > M_SZ) {
      lprintf(stderr, "Unknown component type: (%d) on line %d of file: %s\n",
	      type, *nline, modfile);
      return CMP_ERROR;
    } else {
/*
 * Convert file units to radians.
 */
      radius *= mastor;
      major *= mastor;
      phi *= dtor;
      theta *= dtor;
/*
 * Add the shifted component to the Model.
 */
      if(add_xycmp(mod, docomp, freepar, flux,
		   radius * sin(theta) + east * mastor,
		   radius * cos(theta) + north * mastor,
		   major, ratio, phi, type, freq0, spcind)==NULL)
	return CMP_ERROR;
/*
 * Increment the record of the number of components successfully read.
 */
      return CMP_READ;
    };
  };
  return CMP_EMPTY;
}

/*.......................................................................
 * Write the components of a Model to a file.
 *
 * Input:
 *   mod    Model * The Model to be written. (Can be NULL or empty).
 *   east   float   Eastward X-axis offset to remove before writing (radians).
 *   north  float   Northward Y-axis offset to remove before writing (radians).
 *   docut    int   If not 0 then write only components whose flux
 *                  exceeds the value of 'cut'.
 *   cut    float   See 'docut'.
 *   fd      FILE * The file descriptor of an open file. I am assuming
 *                  that all users of this function will want to start
 *                  the file with a descriptive header before calling
 *                  this routine - hence the requirement for an open
 *                  file rather than the name of a file to be opened.
 * Output:
 *   return   int   0 - No write errors.
 */
int wmodel(Model *mod, float east, float north, int docut, float cut, FILE *fd)
{
  Modcmp *cmp;  /* The component being processed */
/*
 * Label the file columns.
 */
  lprintf(fd, "! Flux (Jy) Radius (mas)  Theta (deg)  Major (mas)  Axial ratio   Phi (deg) T \\\n! Freq (Hz)     SpecIndex\n");
/*
 * No components to be written?
 */
  if(mod==NULL || mod->ncmp==0)
    return 0;
/*
 * Write the components according to their types.
 */
  for(cmp=mod->head; cmp != NULL;cmp=cmp->next) {
/*
 * Apply flux cutoff where requested.
 */
    if(!docut || cmp->flux > cut) {
      double radius; /* Polar radius to component */
      double theta;  /* Polar angle to component */
/*
 * Get the un-shifted X and Y axis positions of the component.
 */
      double xpos = cmp->x - east;
      double ypos = cmp->y - north;
/*
 * Convert component X,Y position to polar representation.
 */
      if(xpos == 0.0 && ypos == 0.0) {
	radius = 0.0;
	theta = 0.0;
      } else {
	radius = rtomas * sqrt(xpos * xpos + ypos * ypos);
	theta  = rtod * atan2(xpos, ypos);
      };
/*
 * All component types require the component flux, radius and theta.
 */
      lprintf(fd,"%#10.6g%c", cmp->flux, cmp->freepar & M_FLUX ?'v':' ');
      lprintf(fd," %#11.6g%c", radius, cmp->freepar & M_CENT ?'v':' ');
      lprintf(fd," %#11.6g%c", theta, cmp->freepar & M_CENT ? 'v':' ');
/*
 * All component types except delta functions also require
 * the major-axis, axis-ratio and axis-direction fields as well.
 * Gaussian components don't need the type field but I will
 * write it anyway for all the above types.
 */
      if(cmp->type != M_DELT || cmp->freq0 > 0.0) {
	lprintf(fd," %#11.6g%c", cmp->major*rtomas,
		cmp->freepar & M_MAJOR ? 'v':' ');
	lprintf(fd," %#11.6g%c", cmp->ratio, cmp->freepar & M_RATIO ? 'v':' ');
	lprintf(fd," %#10.6g%c", cmp->phi*rtod, cmp->freepar & M_PHI ? 'v':' ');
	lprintf(fd," %d", cmp->type);
	if(cmp->freq0 > 0.0) {
	  lprintf(fd, " %s%#11.6g", fd==stdout ? "\\\n ":"", cmp->freq0);
	  lprintf(fd, " %11.6g%c", cmp->spcind,
		  cmp->freepar & M_SPCIND ? 'v':' ');
	};
      };
/*
 * New line.
 */
      lputc('\n', fd);
/*
 * Check for write errors.
 */
      if(ferror(fd)) {
	lprintf(stderr, "Error writing model file\n");
	return -1;
      };
    };
  };
/*
 * Finished successfully.
 */
  return 0;
}

/*.......................................................................
 * Delete all components after the first component that lies below
 * given flux limit. NB This is not the same as removing all components
 * which have a flux below a given level. The model will be deleted if
 * no components remain.
 *
 * Input:
 *   mod  Model *  The model to be modified.
 *   cut  float    The cuttoff level.
 * Output:
 *   return Model * The modified model.
 */
Model *cut_mod(Model *mod, float cut)
{
  Modcmp *cmp;   /* The model component being considered. */
  Modcmp *last;  /* The last un-deleted component in the list */
/*
 * Find the first component that lies below 'cut'. (Keep a record of
 * the 'last' component before the one to be cut so that we know
 * where to terminate the list).
 */
  last=NULL;
  for(cmp=mod->head; cmp != NULL && cmp->flux >= cut; cmp = cmp->next)
    last = cmp;
/*
 * Delete the residual list.
 */
  if(last) {
    while(last->next != NULL)
      del_cmp(rem_cmp(mod,last,last->next));
  };
  return mod;
}

/*.......................................................................
 * Given a model object combine all delta function model components that
 * have the same position, thus shrinking the table.
 *
 * Input/Output:
 *  mod    Model *   The model to be compressed.
 *  return Model *   The squashed model. (same as mod).
 */
Model *squash(Model *mod)
{
  Modcmp *target;    /* The component forming the addition of all similar */
                     /* components */
  Modcmp *cmp;       /* The component being checked against 'target' */
  Modcmp *prev;      /* The component preceding 'cmp' */
/*
 * No-op if no model to be squashed or the model is already squashed.
 */
  if(mod==NULL || mod->issqd)
    return mod;
  mod->issqd = 1;  /* Mark the model as squashed */
/*
 * For each delta-function model component (target) search for any others
 * of the same type with identical positions, add their flux to the target
 * component, unlink them from the component list and delete them.
 * The outer loop gets the next component in the list, then the inner loop
 * searches all components up to that point for matchng positions. If found
 * the outer loop component is unlinked and deleted while the inner loop
 * component has its flux incremented.
 */
  prev = mod->head;
  if(prev!=NULL) {
    for(cmp=prev->next; cmp != NULL; cmp=prev->next) {
      if(cmp->type != M_DELT) {
	prev=prev->next;
	continue;
      };
      for(target=mod->head; target != cmp; target=target->next) {
	if(cmp->x == target->x && cmp->y == target->y && target->type==M_DELT &&
	   cmp->spcind == target->spcind) {
/*
 * Get the reference flux of the component that is to be added.
 */
	  float flux = cmp->flux;
/*
 * If the spectral indeces of the two components (verified as being equal
 * above) are non-zero, get the flux of the component that is being added,
 * at the reference frequency of the target component.
 */
	  if(cmp->spcind != 0.0 && cmp->freq0 != target->freq0)
	    flux *= pow(target->freq0 / cmp->freq0, cmp->spcind);
/*
 * Add the flux of the new component to the target.
 */
	  target->flux += flux;
	  mod->flux += flux;       /* NB. rem_cmp() subtracts removed flux */
/*
 * Unlink and delete the redundant component (cmp).
 */
	  cmp = del_cmp(rem_cmp(mod,prev,cmp));
	  break;
	};
      };
/*
 * If the outer-loop component was not deleted and unlinked then the
 * preceding component in the next iteration is the current component.
 */
      if(target==cmp)
	prev=cmp;
    };
  };
/*
 * Return the squashed model.
 */
  mod->issqd = 1;
  return mod;
}

/*.......................................................................
 * Add one model to another. For example to add 'Model *old' to
 * 'Model *mod', one would type:
 *
 *  mod = add_mod(mod,old,1,1);
 *  old=del_Model(old);
 *
 * If 'old' is NULL then a (squashed if docomp=1) version of 'mod' will be
 * returned directly.
 *
 * Input:
 *  mod   Model *   The model to be appended to. If NULL a new model
 *                  container will be allocated for return, in its place.
 *  old   Model *   The model to be appended to 'mod'. On return this
 *                  model container will be empty and should be deleted
 *                  if no longer required.
 *  docomp int      If true then compress the models together by
 *                  combining delta-function components at identical
 *                  positions.
 *  append int      If true, append the model components of 'old' to
 *                  those of 'mod'. If false (0), prepend the model
 *                  components of 'old' to those of 'mod'.
 * Output:
 *  return Model *  This will be the pointer to 'mod' which will contain
 *                  the combined models, or NULL on error.
 */
Model *add_mod(Model *mod, Model *old, int docomp, int append)
{
  Modcmp *oldcmp; /* A component in 'old' */
  Modcmp *next;   /* The next unprocessed component in 'old' */
/*
 * Allocate an output model container?
 */
  if(mod==NULL) {
    mod = new_Model();
    if(mod==NULL)
      return NULL;
  };
/*
 * If the models are to be compressed together then first compress
 * the individual models.
 */
  if(docomp) {
    squash(mod);
    squash(old);
  };
/*
 * If the 'old' container doesn't exist then simply return the original
 * 'mod' container (squashed if docomp=1).
 */
  if(old==NULL)
    return mod;
/*
 * We only have the tools to append components to an existing model, so
 * if the user requested that the model components of 'old' be prepended
 * to those of 'mod', implement this by first swapping the contents of the
 * two models and then appending the new contents of 'old' to those of
 * 'mod'.
 */
  if(!append) {
    Model tmpmod = *mod;
    *mod = *old;
    *old = tmpmod;
  };
/*
 * Now transfer each model component of 'old' to 'mod'.
 */
  for(oldcmp=old->head; oldcmp != NULL; oldcmp=next) {
    next = oldcmp->next;
    add_cmp(oldcmp, mod, docomp);
  };
/*
 * Mark the old container as empty.
 */
  old->ncmp=0;
  old->issqd = 1;
  old->flux = 0.0;
  old->head=0;
  old->tail=0;
/*
 * Return the combined model.
 */
  return mod;
}

/*.......................................................................
 * Shift centroids of all the components of a model.
 *
 * Input/Output:
 *   mod   Model *  The model to be shifted (can be NULL).
 * Input:
 *   east  float    Shift eastwards  (radians).
 *   north float    Shift northwards (radians).
 */
void shiftmod(Model *mod, float east, float north)
{
  Modcmp *cmp;  /* The model component being shifted */
  if(mod == NULL)
    return;
  for(cmp=mod->head; cmp != NULL; cmp=cmp->next) {
    cmp->x += east;
    cmp->y += north;
  };
  return;
}

/*.......................................................................
 * Append an existing model component to a model.
 *
 * Input:
 *  cmp   Modcmp *  The component to be added.
 *  mod    Model *  The model to be appended to.
 *  docomp   int    If true add the component to any existing component
 *                  with the same X and Y (only delta components). In
 *                  this case 'cmp' itself will be deleted.
 * Output:
 *  return Modcmp * The original 'cmp' or the component to which it was
 *                  merged, or NULL on error.
 */
Modcmp *add_cmp(Modcmp *cmp, Model *mod, int docomp)
{
  Modcmp *oldcmp; /* A component in 'mod' */
/*
 * Check the arguments.
 */
  if(!mod) {
    lprintf(stderr, "NULL model encounterred in add_cmp()\n");
    return del_cmp(cmp);
  };
/*
 * If docmp is true and the component is a delta component, see if
 * there is an existing delta component with the same position
 * to add merge the component with.
 */
  if(docomp && cmp->type==M_DELT) {
    for(oldcmp=mod->head; oldcmp != NULL; oldcmp=oldcmp->next) {
      if(cmp->x == oldcmp->x && cmp->y == oldcmp->y && oldcmp->type==M_DELT &&
	 cmp->spcind == oldcmp->spcind) {
/*
 * Get the reference frequency of the component to be added.
 */
	float flux = cmp->flux;
/*
 * If the equal spectral index of the two components is non-zero,
 * get the flux to be added at the reference frequency of the target
 * component.
 */
	if(cmp->spcind != 0.0 && cmp->freq0 != oldcmp->freq0)
	  flux *= pow(oldcmp->freq0 / cmp->freq0, cmp->spcind);
/*
 * Add the flux to that of the existing component.
 */
	oldcmp->freepar |= cmp->freepar;
	oldcmp->flux += flux;
	mod->flux += flux;
	cmp = del_cmp(cmp);
	return oldcmp;
      };
    };
  };
/*
 * Append the component to the tail of the model component list.
 */
  if(mod->head == NULL) {
    mod->head = mod->tail = cmp;
    mod->issqd = 1;
  }
  else {
    mod->tail->next = cmp;
    mod->tail = cmp;
  };
  mod->ncmp++;
  mod->flux += cmp->flux;
  cmp->next = NULL;
  mod->issqd = mod->issqd && docomp;
/*
 * The first time that a model component is encountered that isn't a
 * continuum delta function, mark the model as containing more than
 * simple delta functions.
 */
  mod->isdelt = mod->isdelt && cmp->type==M_DELT && cmp->spcind==0.0 &&
    !(cmp->freepar & M_SPCIND);
  return cmp;
}

/*.......................................................................
 * Remove a component from a model and return it. In order to unlink the
 * component from the model component list, the function needs to 
 * know what the previous component in the list is. If you don't know
 * this set prev=NULL and a linear search will be made through the list
 * for the previous component. Needless to say, if you are iterating
 * over the model list it is faster to keep track of the previous
 * unremoved component than to get this function to search for it.
 * However, you must set prev=NULL for the first component since there
 * is no preceding component in the list.
 *
 * Input:
 *  mod     Model *  The model containing the component to be removed.
 *  prev   Modcmp *  The component in the model list preceding 'cmp', or
 *                   NULL if unknown or cmp is the head of the list.
 *  cmp    Modcmp *  The model component to be removed.
 * Output:
 *  return Modcmp *  The extracted component. Apply del_cmp() to this
 *                   if it is no longer required.
 */
Modcmp *rem_cmp(Model *mod, Modcmp *prev, Modcmp *cmp)
{
  static Modcmp *tmpcmp;
  if(cmp==NULL)
    return cmp;
/*
 * If prev is NULL search for 'cmp' and the component preceding it.
 */
  if(prev==NULL) {
    for(tmpcmp=mod->head; tmpcmp!=cmp && tmpcmp!=NULL; tmpcmp=tmpcmp->next)
      prev = tmpcmp;
    if(tmpcmp==NULL) {
      lprintf(stderr, "rem_cmp: programmer error: component not found\n");
      return cmp;
    };
  };
/*
 * If 'prev' is still NULL then 'cmp' is the first component of the model.
 */
  if(prev==NULL)
    mod->head = cmp->next;
  else
    prev->next = cmp->next;
/*
 * Fix up the tail if cmp was the tail of the list.
 */
  if(cmp->next==NULL)
    mod->tail = prev;
/*
 * Record the shortening of the list.
 */
  mod->ncmp--;
  mod->flux -= cmp->flux;
/*
 * Reset the model characteristics if the model is now empty.
 */
  if(mod->ncmp==0) {
    mod->issqd = mod->isdelt = 1;
    mod->flux = 0;
  };
  return cmp;
}

/*.......................................................................
 * Re-arrange the contents of two models such that one contains just
 * fixed components and the other contains just variable components. The
 * later is intended for use in model fitting.
 *
 * Input:
 *  amod    Model * The Model to be returned containing the fixed part
 *                  of the model. This must not be NULL.
 *  bmod    Model * The Model to be returned containing the variable part
 *                  of the model. This must not be NULL.
 * Output:
 *  return    int   0 - OK.
 *                  1 - Error.
 */
int var_mod(Model *amod, Model *bmod)
{
  Modcmp *cmp;   /* The model component being considered. */
  Modcmp *next;  /* The next model component in the list after 'cmp' */
  Modcmp *last;  /* The last unremoved component in the list before 'cmp' */
/*
 * We require two model containers.
 */
  if(amod==NULL || bmod==NULL) {
    lprintf(stderr, "var_mod: NULL Model container intercepted.\n");
    return 1;
  };
/*
 * Loop through the components of 'amod', remove variable components
 * and append them to 'bmod'. 'last' is the last unremoved component
 * before the current one, 'cmp' and 'next' is the next one after it in
 * the list. 
 */
  last = NULL;
  for(cmp=amod->head; cmp != NULL; cmp=next) {
    next = cmp->next;
    if(!cmp->freepar)
      last=cmp;
    else      /* Move variable components to 'bmod' */
      cmp=add_cmp(rem_cmp(amod,last,cmp), bmod, 1);
  };
/*
 * Loop through the components of 'bmod', remove fixed components
 * and append them to 'amod'. 'last' is the last unremoved component
 * before the current one, 'cmp' and 'next' is the next one after it in
 * the list. 
 */
  last = NULL;
  for(cmp=bmod->head; cmp != NULL; cmp=next) {
    next = cmp->next;
    if(cmp->freepar)
      last=cmp;
    else      /* Move fixed components to 'amod' */
      cmp=add_cmp(rem_cmp(bmod,last,cmp), amod, 1);
  };
  return 0;
}

static int ed_wmod(const char *modfile, Model *mod);

/*.......................................................................
 * Alow the user to edit a model in an external editor.
 *
 * Input:
 *  mod     Model *   The model to be edited (NULL is legal).
 * Output:
 *  return  Model *   The edited model, in which case 'mod' will have been
 *                    deleted, or the undeleted 'mod' on error.
 *                    Use like:  model = ed_model(model);
 */
Model *ed_model(Model *mod)
{
  char *modfile;  /* The name of the model scratch file */
  int edited=0;   /* True when the model has been sucessfully edited */
/*
 * Get the name of a scratch file within which to write the current model.
 */
  modfile = scrname("edmod.scr");
/*
 * Write the existing contents of the model to the scratch file,
 * and allow the user to edit the file.
 */
  if(modfile) {
    if(ed_wmod(modfile, mod)==0 && ed_file(modfile)==0) {
/*
 * Create a new model container to read the modified model into.
 */
      Model *newmod = new_Model();
/*
 * If a new container was sucessfully allocated, fill it from the
 * edited model file, and replace the old model with it.
 */
      if(newmod==NULL || rmodel(newmod, 0.0f, 0.0f, 1, modfile)) {
	newmod = del_Model(newmod);
      } else {
	del_Model(mod);
	mod = newmod;
	edited = 1;   /* Mark the model as sucesfully edited */
      };
    };
/*
 * Delete the model scratch file.
 */
    remove(modfile);
/*
 * Release the memory used to store the scratch file name.
 */
    free(modfile);
  };
/*
 * Return the original or edited model.
 */
  if(!edited)
    lprintf(stdout, "Reinstating the original un-edited model.\n");
  return mod;
}

/*.......................................................................
 * Private function of ed_model(), used to write the components of a
 * model to a given scratch file.
 *
 * Input:
 *  modfile     char *   The name for the model scratch file.
 *  mod        Model *   The model to be written.
 * Output:
 *  return       int     0 - OK.
 *                       1 - Error.
 */
static int ed_wmod(const char *modfile, Model *mod)
{
/*
 * Create and open the scratch file for writing.
 */
  FILE *fp = fopen(modfile, "w");
  if(fp == NULL) {
    lprintf(stderr, "ed_model: Unable to open scratch file: %s\n", modfile);
    return 1;
  };
/*
 * Write the contents of the model to the scratch file.
 */
  wmodel(mod, 0.0f, 0.0f, 0, 0.0f, fp);
/*
 * Close the model file.
 */
  if(fclose(fp)==EOF)
    return 1;
  return 0;
}

/*.......................................................................
 * Allocate a new copy of a given model.
 *
 * Input:
 *  mod       Model *   The model to be copied.
 * Output:
 *  return    Model *   The newly allocated copy of the model, or NULL
 *                      on error.
 */
Model *cpy_Model(Model *mod)
{
  Model *newmod;   /* The new model */
  Modcmp *cmp;     /* A model component being copied. */
/*
 * Copy a NULL model?
 */
  if(!mod)
    return NULL;
/*
 * Create the new model container.
 */
  newmod = new_Model();
  if(!newmod)
    return NULL;
/*
 * Copy the old model members.
 */
  *newmod = *mod;
/*
 * Initialize the parts that should be assigned anew.
 */
  newmod->head = NULL;
  newmod->tail = NULL;
/*
 * Add copies of the original components to the new model.
 */
  for(cmp=mod->head; cmp; cmp=cmp->next) {
/*
 * Allocate a new component.
 */
    Modcmp *newcmp = new_cmp(++newmod->ncmp);
    if(!newcmp)
      return del_Model(newmod);
/*
 * Copy the old component into the new component.
 */
    *newcmp = *cmp;
/*
 * Append the new component.
 */
    newcmp->next = NULL;
    if(newmod->tail) {
      newmod->tail->next = newcmp;
      newmod->tail = newcmp;
    } else {
      newmod->head = newmod->tail = newcmp;
    };
  };
  return newmod;
}
