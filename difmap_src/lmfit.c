#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matinv.h"
#include "lmfit.h"

static Lmfit *lm_nomem(Lmfit *lm);
static int lm_getfit(Lmfit *lm);

static void clr_Fitpar(Fitpar *fp);
static void cop_Fitpar(Fitpar *to, Fitpar *from, int nfree);
static Fitpar *new_Fitpar(Fitpar *fp, int nfree);
static Fitpar *del_Fitpar(Fitpar *fp);

/*.......................................................................
 * Create and initialize a Levenberg-Marquardt fit object.
 *
 * Input:
 *  obj     void *  The pointer to an anonymous object used by the
 *                  user functions given below, to access
 *                  the model and data being fitted.
 *  nfree    int    The number of free parameters to be fitted.
 *  getfree         A function (declared as GETFREE(getfree)), whos
 *                  function is to return the current values of the model
 *                  free parameters in a provided array.
 *                    Input arguments:
 *                      obj    void *  The value of obj above will be sent.
 *                      nfree   int    The value of nfree above will be sent.
 *                    Input/Output argument:
 *                      pars double *  An array of nfree elements will be
 *                                     provided. On return this should
 *                                     contain the values of each free
 *                                     parameter in the model.
 *                    Output argument:
 *                      return  int    0 - OK.
 *                                     1 - Error.
 *
 *  setfree         A function (declared as SETFREE(setfree)), whos
 *                  function is to assign given values for the model free
 *                  parameters to the model.
 *                    Input arguments:
 *                      obj    void *  The value of obj above will be sent.
 *                      nfree   int    The value of nfree above will be sent.
 *                      pars double *  An array of nfree elements will be
 *                                     provided cotaining the new values for
 *                                     the model free parameters. These should
 *                                     be recorded in the model.
 *                    Output argument:
 *                      return  int    0 - OK.
 *                                     1 - Error.
 *
 *  getnext         A function (declared as GETNEXT(getnext)), who's
 *                  function is to return information about a single
 *                  data point. Each call must return details about a
 *                  different data point until all data points have
 *                  been seen. On the following call 0 should be
 *                  returned to signify end-of-data, and on call
 *                  following that call, the cycle should restart with
 *                  the first data point. The whole iteration cycle
 *                  will be invoked each time lm_fit() is called.
 *                    Input arguments:
 *                      obj     void * The value of obj above will be sent.
 *                      nfree    int   The value of nfree above will be sent.
 *                    Input/Output arguments:
 *                      dy    double * *dy must be assigned with the
 *                                     difference between the model
 *                                     and the data, as (data - model).
 *                      wt   double * *wt must be assigned the weight of
 *                                     the measurement, defined as the
 *                                     reciprocal of the square of the
 *                                     predicted standard deviation.
 *                      mgrad double * An array of nfree elements will
 *                                     be provided. Each element must
 *                                     be assigned with the partial
 *                                     derivative of the model wrt the
 *                                     corresponding free parameter.
 *                    Output arguments:
 *                      return   int   0 - End of data.
 *                                     1 - Valid data point returned.
 *                                    -1 - Error.
 * Output:
 *  return  Lmfit * The initialized fit object, or NULL on error.
 */
Lmfit *new_Lmfit(void *obj, int nfree, GETFREE(*getfree), SETFREE(*setfree),
                 GETNEXT(*getnext))
{
  Lmfit *lm;   /* object to be returned */
/*
 * Check arguments.
 */
  if(nfree<=0) {
    fprintf(stderr, "new_Lmfit: nfree must be > 0.\n");
    return NULL;
  };
  if(getfree==0 || setfree==0 || getnext==0) {
    fprintf(stderr, "new_Lmfit: One or more access functions given as NULL.\n");
    return NULL;
  };
/*
 * Allocate the container.
 */
  lm = (Lmfit *) malloc(sizeof(Lmfit));
  if(lm==NULL)
    return lm_nomem(lm);
/*
 * Initialize it at least to the point at which it can safely be sent to
 * del_Lmfit().
 */
  lm->nfree = nfree;
  lm->obj = obj;
  clr_Fitpar(&lm->best);
  clr_Fitpar(&lm->new);
  lm->work = NULL;
  lm->dpwork = NULL;
  lm->iwork1 = NULL;
  lm->iwork2 = NULL;
  lm->iwork3 = NULL;
  lm->incfac = 0.001;
  lm->getfree = getfree;
  lm->setfree = setfree;
  lm->getnext = getnext;
/*
 * Allocate the contents of the fit containers.
 */
  if(new_Fitpar(&lm->best, nfree)==NULL || new_Fitpar(&lm->new, nfree)==NULL)
    return lm_nomem(lm);
/*
 * Allocate 1D work arrays.
 */
  lm->work = (double *) malloc(sizeof(double) * nfree);
  lm->dpwork = (double **) malloc(sizeof(double *) * nfree);
  lm->iwork1 = (int *) malloc(sizeof(int) * nfree);
  lm->iwork2 = (int *) malloc(sizeof(int) * nfree);
  lm->iwork3 = (int *) malloc(sizeof(int) * nfree);
  if(lm->work==NULL || lm->dpwork==NULL || lm->iwork1==NULL ||
     lm->iwork2==NULL || lm->iwork3==NULL)
    return lm_nomem(lm);
/*
 * Load the initial values of the free parameters into the
 * array of trial parameter values.
 */
  if(getfree(lm->obj, lm->nfree, lm->new.pars))
    return del_Lmfit(lm);
/*
 * Assign an impossible value as the best-fit chi-squared value, to
 * signal that there is no best fit yet.
 */
  lm->best.chisq = -1.0;
/*
 * Return the initialized object.
 */
  return lm;
}

/*.......................................................................
 * Private error return function of new_Lmfit().
 */
static Lmfit *lm_nomem(Lmfit *lm)
{
  fprintf(stderr, "new_Lmfit: Insufficient memory.\n");
  return del_Lmfit(lm);
}

/*.......................................................................
 * Clear the un-allocated contents of a container of fit parameters,
 * such that del_Fitpar() can detect that nothing has been allocated.
 *
 * Input:
 *  fp    Fitpar *  The container to be cleared.
 */
static void clr_Fitpar(Fitpar *fp)
{
  fp->hessian = NULL;
  fp->cgrad = NULL;
  fp->pars = NULL;
  fp->chisq = 0.0;
  fp->rchisq = 0.0;
  fp->ndfree = 0;
}

/*.......................................................................
 * Allocate the contents of a container of fit parameters.
 *
 * Input:
 *  fp     Fitpar *  The container to be filled.
 *  nfree     int    The number of free parameters in the model.
 * Output:
 *  return Fitpar *  The same as 'fp', or NULL on error.
 */
static Fitpar *new_Fitpar(Fitpar *fp, int nfree)
{
  int row;   /* Hessian matrix row index */
/*
 * Allocate the hessian matrix dope vector.
 */
  fp->hessian = (double **) malloc(sizeof(double *) * nfree);
  if(fp->hessian==NULL)
    return NULL;
/*
 * Allocate the body of the matrix.
 */
  fp->hessian[0] = (double *) malloc(sizeof(double) * nfree * nfree);
  if(fp->hessian[0]==NULL)
    return NULL;
/*
 * Initialize the dope vector thread the matrix.
 */
  for(row=1; row<nfree; row++)
    fp->hessian[row] = fp->hessian[row-1] + nfree;
/*
 * Allocate the arrays used to record the Chi-squared gradient and parameter
 * values.
 */
  fp->cgrad = (double *) malloc(sizeof(double) * nfree);
  fp->pars = (double *) malloc(sizeof(double) * nfree);
  if(fp->cgrad==NULL || fp->pars == NULL)
    return NULL;
/*
 * Return the filled container.
 */
  return fp;
}

/*.......................................................................
 * Copy the contents of one container of fit parameters to another.
 *
 * Input:
 *  dest   Fitpar *  The container to copy into.
 *  from   Fitpar *  The container to copy from.
 *  nfree     int    The number of free parameters in the model.
 */
static void cop_Fitpar(Fitpar *dest, Fitpar *from, int nfree)
{
  int row;  /* The row being copied */
  int col;  /* The column being copied */
/*
 * Copy the hessian matrix.
 */
  for(row=0; row<nfree; row++) {
    for(col=0; col<nfree; col++) {
      dest->hessian[row][col] = from->hessian[row][col];
    };
  };
/*
 * Copy the Chi-squared-gradient and free-parameter arrays.
 */
  for(col=0; col<nfree; col++) {
    dest->cgrad[col] = from->cgrad[col];
    dest->pars[col] = from->pars[col];
  };
/*
 * Copy chi-squared.
 */
  dest->chisq = from->chisq;
  dest->rchisq = from->rchisq;
  dest->ndfree = from->ndfree;
  return;
}

/*.......................................................................
 * Delete the contents of a container of fit parameters.
 *
 * Input:
 *  fp     Fitpar *  The container to be emptied.
 * Output:
 *  return Fitpar *  The deleted container, ie. NULL.
 */
static Fitpar *del_Fitpar(Fitpar *fp)
{
  if(fp) {
    if(fp->hessian) {
      if(*fp->hessian)
	free(*fp->hessian);  /* Body of matrix */
      free(fp->hessian);     /* Dope vector */
    };
    if(fp->cgrad)
      free(fp->cgrad);
    if(fp->pars)
      free(fp->pars);
  };
  return NULL;
}

/*.......................................................................
 * Delete a LM fit object.
 *
 * Input:
 *  lm      Lmfit *  The instance to be deleted.
 * Output:
 *  return  Lmfit *  The deleted object, ie. NULL.
 *                   Use like: lm = del_Lmfit(lm);
 */
Lmfit *del_Lmfit(Lmfit *lm)
{
  if(lm) {
    del_Fitpar(&lm->best);
    del_Fitpar(&lm->new);
    if(lm->work)
      free(lm->work);
    if(lm->dpwork)
      free(lm->dpwork);
    if(lm->iwork1)
      free(lm->iwork1);
    if(lm->iwork2)
      free(lm->iwork2);
    if(lm->iwork3)
      free(lm->iwork3);
/*
 * Free the descriptor container.
 */
    free(lm);
  };
  return NULL;
}

/*.......................................................................
 * Perform a single iteration of the Levenberg-Marquardt minimization
 * non-linear least-squares technique.
 *
 * On return the parameters that gave the best fit in the current or
 * previous iterations will be re-established in the input model.
 *
 * Input:
 *  lm       Lmfit *  A fit descriptor initially returned by new_Lmfit().
 * Output:
 *  return Lmstate    LM_ABORT  - Error - no further iteration is possible.
 *                                An error message will have been presented
 *                                to the user, explaining why.
 *                    LM_BETTER - An improved fit is described in lm->best.
 *                    LM_WORSE  - The trial fit (now described in lm->new)
 *                                failed to improve over the previous best
 *                                fit (still described in lm->best). This
 *                                does not mean that you should stop iterating.
 *                                The next iteration will re-try with more
 *                                conservative parameter increments.
 */
Lmstate lm_fit(Lmfit *lm)
{
  Lmstate state;  /* The return fit-status code */
/*
 * Get the new chi-squared, the new hessian matrix, and model vs. free
 * parameter gradients for the trial parameter values in lm->new.pars.
 * The values are returned in lm->new{}.
 */
  if(lm_getfit(lm)) {
    state = LM_ABORT;        /* Error */
/*
 * If the reduced chi-squared decreased or this is the first
 * iteration (marked by chisq<=0.0), record the new fit as the best fit.
 */
  } else if(lm->new.rchisq < lm->best.rchisq || lm->best.chisq <= 0.0) {
    state = LM_BETTER;
    cop_Fitpar(&lm->best, &lm->new, lm->nfree);
/*
 * Try a coarser parameter increment on the next iteration.
 */
    lm->incfac *= 0.5;
/*
 * We stepped too far so try a finer parameter increment on the next
 * iteration.
 */
  } else {
    state = LM_WORSE;
    lm->incfac *= 10.0;
  };
/*
 * Ensure that the user's model always contains the best-fit parameters
 * at the end of each iteration.
 */
  if(lm->setfree(lm->obj, lm->nfree, lm->best.pars))
    state = LM_ABORT;
/*
 * Return an indication of whether the trial fit succeded.
 */
  return state;
}

/*.......................................................................
 * Iterate through the data vs. model in order to get the hessian matrix,
 * Chi-squared vs. free-parameter gradients and the current value of
 * chi-squared.
 * The fit details will be returned in lm->new{}.
 *
 * Input:
 *  lm       Lmfit *   The fit descriptor.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
static int lm_getfit(Lmfit *lm)
{
  int iret;     /* Return code from user getnext() function */
  int row;      /* A row index */
  int col;      /* A column index */
  double dy;    /* Difference between model and data as (data-model) */
  double wt;    /* The data weight as 1/(standard_deviation^2). */
/*
 * If this is not the first iteration, use the Hessian matrix from the
 * last best fit to solve for the new trial parameters.
 */
  if(lm->best.chisq > 0.0) {   /* Not the first iteration? */
/*
 * Get the best-fit hessian.
 */
    cop_Fitpar(&lm->new, &lm->best, lm->nfree);
/*
 * Exagerate the step increments by introducing the extra incfac scale
 * factor on the diagonal of the Hessian matrix.
 */
    for(col=0; col<lm->nfree; col++)
      lm->new.hessian[col][col] *= (1.0 + lm->incfac);
/*
 * Solve for the free-parameter increments via LU decomposition
 * and back-substitution.
 */
    if(gj_solve(lm->new.hessian,lm->new.cgrad,lm->iwork1,lm->dpwork,lm->nfree))
      return 0;
/*
 * Apply the parameter increments (now in lm->new.cgrad) to prepare
 * trial parameters for the next iteration. 
 */
    for(col=0; col<lm->nfree; col++)
      lm->new.pars[col] = lm->best.pars[col] + lm->new.cgrad[col];
/*
 * Establish the new trial parameters in the user's model.
 */
    if(lm->setfree(lm->obj, lm->nfree, lm->new.pars))
      return 1;
  };
/*
 * Get the trial free parameter values. Even though setfree() may have been
 * called above, this is necessary since setfree() may have limited the
 * values set.
 */
  if(lm->getfree(lm->obj, lm->nfree, lm->new.pars))
    return 1;
/*
 * Initialize pertinent members of the output fit descriptor.
 */
  for(row=0; row<lm->nfree; row++) {
    lm->new.cgrad[row] = 0.0;
    for(col=0; col<lm->nfree; col++) {
      lm->new.hessian[row][col] = 0.0;
    };
  };
/*
 * The number of degrees of freedom is given by the number of measurements
 * minus the number of free parameters. We know the number of free parameters
 * now so initialize with this. We will then add 1 to this parameter as
 * each measurement is processed, and if then lm->ndfree has not become
 * +ve then we know that the model is too complex.
 */
  lm->new.ndfree = -lm->nfree;
/*
 * Initialize the chi-squared sum.
 */
  lm->new.chisq = 0.0;
/*
 * Sum over data and model.
 */
  while((iret=lm->getnext(lm->obj, lm->nfree, &dy, &wt, lm->work))==1) {
    double *mgrad = lm->work; /* Model gradient */
/*
 * Record the receipt of a new measurement.
 */
    lm->new.ndfree++;
/*
 * Add in the contribution to the Hessian matrix, Chi-squared gradient
 * and chi-squared, of the latest measurement.
 */
    for(row=0; row<lm->nfree; row++) {
      double tmp = wt * mgrad[row];
      double *hessian_row = lm->new.hessian[row];
      for(col=0; col<=row; col++) {
	hessian_row[col] += tmp * mgrad[col];
      };
      lm->new.cgrad[row] += dy * tmp;
    };
    lm->new.chisq += wt * dy * dy;
  };
/*
 * Fill in the symmetric upper half of the Hessian matrix.
 */
  if(iret==0) {
    for(row=0; row<lm->nfree; row++) {
      for(col=row+1; col<lm->nfree; col++) {
	lm->new.hessian[row][col] = lm->new.hessian[col][row];
      };
    };
/*
 * Determine the value of reduced chi-squared.
 */
    if(lm->new.ndfree >= 1) {
      lm->new.rchisq = lm->new.chisq / lm->new.ndfree;
    } else {
      fprintf(stderr, "lm_getfit: Fewer measurements than free parameters.\n");
      iret = 1;
    };
  };
  return iret != 0;
}

/*.......................................................................
 * Return the covariance matrix of the current fit.
 *
 * Input:
 *  lm      Lmfit *  The fit descriptor.
 * Output:
 *  return double ** The covariance matrix, or NULL on error.
 *                   The matrix is organised as an array of lm->nfree
 *                   pointers to row arrays of lm->nfree columns. The
 *                   returned matrix is actually work matrix in lm used
 *                   for other purposes, and so may be changed on the
 *                   next call to any other lm_*() function, and will be
 *                   deleted when del_Lmfit(lm) is called.
 */
double **lm_covar(Lmfit *lm)
{
/*
 * Copy the best-fit hessian to lm->new.
 */
  cop_Fitpar(&lm->new, &lm->best, lm->nfree);
/*
 * Invert the copy of the hessian matrix to produce the covariance
 * matrix.
 */
  if(gj_invert(lm->new.hessian, lm->iwork1, lm->iwork2, lm->iwork3, lm->nfree))
    return NULL;
/*
 * Return the covariance matrix.
 */
  return lm->new.hessian;
}

