#ifndef lmfit_h
#define lmfit_h

#define GETFREE(fn) int (fn)(void *obj, int nfree, double *pars)
#define SETFREE(fn) int (fn)(void *obj, int nfree, double *pars)
#define GETNEXT(fn) int (fn)(void *obj, int nfree, double *dy, double *wt, double *mgrad)

typedef struct {
  double **hessian; /* Linearized hessian matrix: (*hessian[nfree])[nfree] */
  double *cgrad;    /* Chi-squared gradient wrt each free parameter */
  double *pars;     /* Array of the nfree latest trial parameter values */
  double chisq;     /* Chi-squared value attained */
  double rchisq;    /* Reduced chi-squared = (chisq / ndfree) */
  long ndfree;      /* The number of degrees of freedom, given as the */
                    /* number of measurements minus lm->nfree */
} Fitpar;

typedef struct {
  int nfree;        /* The number of free parameters */
  void *obj;        /* Pointer to the data/iterator container object */
  Fitpar best;      /* The details of the best fit attained so far */
  Fitpar new;       /* The details of the latest trial fit */
  double *work;     /* Temporary work array of nfree elements */
  double **dpwork;  /* Work dope vector of nfree elements */
  int *iwork1;      /* Temporary work array of nfree elements */
  int *iwork2;      /* Temporary work array of nfree elements */
  int *iwork3;      /* Temporary work array of nfree elements */
  double incfac;    /* The dynamic factor magnifying step sizes */
  GETFREE(*getfree);/* Function to get a copy of the current free parameters */
  SETFREE(*setfree);/* Function to set the current model free parameters */
  GETNEXT(*getnext);/* Get the next data - model,derivatives from iterator */
} Lmfit;

/* Construct and intialize a new Levenberg-Marquardt fit object */

Lmfit *new_Lmfit(void *obj, int nfree, GETFREE(*getfree), SETFREE(*setfree),
		 GETNEXT(*getnext));

/* Delete a Levenberg-Marquardt fit object */

Lmfit *del_Lmfit(Lmfit *lm);

/* Perform a single iteration of minimization - return chi-squared */

typedef enum {LM_ABORT, LM_BETTER, LM_WORSE} Lmstate;

Lmstate lm_fit(Lmfit *lm);

/* Return the covariance matrix of the current best fit */

double **lm_covar(Lmfit *lm);

#endif
