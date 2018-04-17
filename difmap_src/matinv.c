#include <stdio.h>
#include <math.h>

#include "matinv.h"

#define TINY 1.0e-20

/*.......................................................................
 * Use Gauss-Jordan elimination to invert a matrix.
 *
 * Input:
 *  lhs  double **  A SQUARE matrix organised as a C array of pointers
 *                  to vectors, and dimensioned as lhs[nu][nu].
 *  work1   int *   Work array of length nu.
 *  work2   int *   Work array of length nu.
 *  work3   int *   Work array of length nu.
 *  nu      int     The number of unknowns and thus the dimensions of
 *                  lhs and rhs above.
 * Output:
 *  return int      0 - success, anything else means an error occurred.
 *                 -1 - singular matrix.
 */
int gj_invert(double **lhs, int *work1, int *work2, int *work3, int nu)
{
  int *done=work1;   /* Records which diagonals have been pivotted */
  int *orig=work2;   /* Records original row index of n'th pivot row */
  int *dest=work3;   /* Records destination row index of n'th pivot row */
  int iter;          /* Iteration counter. Every iteration reduces one row */
  int maxcol=0;      /* The diagonal element, lhs[diag][diag] being made 1 */
  int maxrow=0;      /* Row with new pivot in it, before row swap. */
  int icol;          /* Column number */
  int irow;          /* Row number */
  double maxpiv=0.0; /* Maximum (absolute) usable pivot element in 'lhs'. */
  double pivval=0.0; /* Signed value of maxpiv */
/*
 * No diagonals yet pivotted.
 */
  for(irow=0; irow<nu; irow++)
    done[irow]=0;
/*
 * Now loop to reduce all diagonal elements to 1 and all other elements
 * to zero. NB the rows and columns and done in the order for which
 * pivot elements are found.
 */
  for(iter=0; iter < nu; iter++) {
/*
 * Find the largest element of the array for a diagonal pivot element
 * that hasn't been pivotted yet.
 */
    maxpiv = 0.0;  /* Records max usable (absolute) value in 'lhs' */
    for(irow=0; irow<nu; irow++) {
      if(!done[irow]) {
	for(icol=0; icol<nu; icol++) {
	  if(!done[icol]) {
	    double newval = lhs[irow][icol];
	    double absval = newval >= 0.0 ? newval : -newval;
	    if(absval > maxpiv) {
	      pivval = newval;
	      maxpiv = absval;
	      maxrow = irow;
	      maxcol = icol;
	    };
	  };
	};
      };
    };
/*
 * If no appropriate pivot was found, the matrix is singular
 * and so can't be solved.
 */
    if(maxpiv==0.0) {
      fprintf(stderr, "gj_invert: Singular matrix.\n");
      return -1;
    };
/*
 * Flag the diagonal that the newly found pivot will be used on, as done.
 */
    done[maxcol] = 1;
/*
 * Swap rows so that the new pivot will be on the associated diagonal.
 */
    if(maxrow != maxcol) {
      double *row1 = lhs[maxrow];
      double *row2 = lhs[maxcol];
      for(icol=0; icol<nu; icol++) {
	double dtmp = row1[icol];
	row1[icol] = row2[icol];
	row2[icol] = dtmp;
      };
    };
/*
 * The row swap results in a column swap in the inverse matrix.
 * In order to be able to unscamble the final result, record the
 * origin and destination of the pivot row here.
 */
    orig[iter] = maxrow;
    dest[iter] = maxcol;
/*
 * Divide the pivot row by its pivot element.
 * The pivot element is redundant after this point, and will not be required
 * hereafter, so before the division, substitute the 1 of the identity matrix
 * such that the associated element of the inverse matrix will be built up.
 */
    lhs[maxcol][maxcol] = 1.0;
    for(icol=0; icol<nu; icol++)
      lhs[maxcol][icol] /= pivval;
/*
 * Now subtract the appropriate amount from all other rows
 * to make them 0 on the column of the pivot element.
 * Start with the rows above the pivot element.
 */
    for(irow=0; irow<maxcol; irow++) {
      double dtmp = lhs[irow][maxcol];
      lhs[irow][maxcol] = 0.0;        /* Substitute identity matrix */
      for(icol=0; icol < nu; icol++)
	lhs[irow][icol] -= lhs[maxcol][icol]*dtmp;
    };
/*
 * Reduce the rows above the pivot row.
 */
    for(irow=maxcol+1; irow<nu; irow++) {
      double dtmp = lhs[irow][maxcol];
      lhs[irow][maxcol] = 0.0;        /* Substitute identity matrix */
      for(icol=0; icol < nu; icol++)
	lhs[irow][icol] -= lhs[maxcol][icol]*dtmp;
    };
  };
/*
 * Now unscramble the column swaps in the inverse matrix that were engendered
 * by pivot row swaps.
 */
  for(iter=nu-1; iter>=0; iter--) {
    int acol = dest[iter];
    int bcol = orig[iter];
    for(irow=0; irow<nu; irow++) {
      double dtmp = lhs[irow][acol];
      lhs[irow][acol] = lhs[irow][bcol];
      lhs[irow][bcol] = dtmp;
    };
  };
  return 0;
}

/*.......................................................................
 * Use Gauss-Jordan elimination to solve a set of linear equations
 * without the overhead of recording the implictly inverted unit matrix.
 *
 * Input:
 *  lhs   double **  A SQUARE matrix organised as a C array of pointers
 *                   to vectors, and dimensioned as lhs[nu][nu].
 *                   NB. The contents of the rows of the matrix will be
 *                       returned in a garbled form.
 *  rhs   double *   The vector of nu right-hand-side values of the
 *                   nu simultaneous equations being solved. On output
 *                   this will contain the nu solutions.
 *  work1    int *   Work array of nu (int)'s.
 *  work2 double **  Work array of nu pointers to (double *).
 *  nu       int     The number of unknowns and thus the dimensions of
 *                   lhs and rhs above.
 * Output:
 *  return   int      0 - success, anything else means an error occurred.
 *                   -1 - singular matrix.
 */
int gj_solve(double **lhs, double *rhs, int *work1, double **work2, int nu)
{
  int *done=work1;   /* Records which diagonals have been pivotted */
  double **rows = work2; /* A row permuted copy of rhs[] */
  int iter;          /* Iteration counter. Every iteration reduces one row */
  int maxcol=0;      /* The diagonal element, lhs[diag][diag] being made 1 */
  int maxrow=0;      /* Row with new pivot in it, before row swap. */
  int icol;          /* Column number */
  int irow;          /* Row number */
  double maxpiv=0.0; /* Maximum (absolute) usable pivot element in 'rows'. */
  double pivval=0.0; /* Signed value of maxpiv */
/*
 * Copy the dope vector of lhs[] and mark all rows as unpivotted.
 */
  for(irow=0; irow<nu; irow++) {
    rows[irow] = lhs[irow];
    done[irow] = 0;
  };
/*
 * Now loop to reduce all diagonal elements to 1 and all other elements
 * to zero. NB the rows and columns and done in the order for which
 * pivot elements are found.
 */
  for(iter=0; iter < nu; iter++) {
/*
 * Find the largest element of the array for a diagonal pivot element
 * that hasn't been pivotted yet.
 */
    maxpiv = 0.0;  /* Records max usable (absolute) value in 'rows' */
    for(irow=0; irow<nu; irow++) {
      if(!done[irow]) {
	for(icol=0; icol<nu; icol++) {
	  if(!done[icol]) {
	    double newval = rows[irow][icol];
	    double absval = newval >= 0.0 ? newval : -newval;
	    if(absval > maxpiv) {
	      pivval = newval;
	      maxpiv = absval;
	      maxrow = irow;
	      maxcol = icol;
	    };
	  };
	};
      };
    };
/*
 * If no appropriate pivot was found, the matrix is singular
 * and so can't be solved.
 */
    if(maxpiv==0.0) {
      fprintf(stderr, "gj_solve: Singular matrix.\n");
      return -1;
    };
/*
 * Flag the diagonal that the newly found pivot will be used on, as done.
 */
    done[maxcol] = 1;
/*
 * Swap rows so that the new pivot will be on the associated diagonal.
 * Simply swap row pointers in the row[] copy of rhs[].
 */
    if(maxrow != maxcol) {
      {                              /* Swap lhs rows */
	double *tmp = rows[maxrow];
	rows[maxrow] = rows[maxcol];
	rows[maxcol] = tmp;
      };
      {                              /* Swap rhs rows */
	double tmp = rhs[maxrow];
	rhs[maxrow] = rhs[maxcol];
	rhs[maxcol] = tmp;
      };
    };
/*
 * Divide the pivot row by its pivot element.
 */
    for(icol=0; icol<nu; icol++)
      rows[maxcol][icol] /= pivval;
    rhs[maxcol] /= pivval;
/*
 * Now subtract the appropriate amount from all other rows
 * to make them 0 on the column of the pivot element.
 * Start with the rows above the pivot element.
 */
    for(irow=0; irow<maxcol; irow++) {
      double dtmp = rows[irow][maxcol];
      for(icol=0; icol < nu; icol++)
	rows[irow][icol] -= rows[maxcol][icol]*dtmp;
      rhs[irow] -= rhs[maxcol]*dtmp;
    };
/*
 * Reduce the rows above the pivot row.
 */
    for(irow=maxcol+1; irow<nu; irow++) {
      double dtmp = rows[irow][maxcol];
      for(icol=0; icol < nu; icol++)
	rows[irow][icol] -= rows[maxcol][icol]*dtmp;
      rhs[irow] -= rhs[maxcol]*dtmp;
    };
  };
  return 0;
}

/*.......................................................................
 * Solve a set of linear simultaneous equations via LU decomposition
 * and backsubstitution, with the option of iterative improvement.
 *
 * Input:
 *  lhs  double **  A SQUARE matrix organised as a C array of pointers
 *                  to vectors, and dimensioned as lhs[nu][nu].
 *                  This is the matrix of known A coefficients on the left
 *                  hand side of the matrix equation  A_ij.x_j = b_j, where
 *                  x_j is the jth unknown and b_j is the j'th known
 *                  coefficient that is independant of any x_j.
 *  rhs  double *   The right-hand side array of the equation detailed
 *                  in the above documentation of 'lhs'. It has dimension
 *                  rhs[nu].
 *  iwrk    int *   Work array of length nu. This could have been allocated
 *                  within this function, but the redundant malloc() free()
 *                  overhead for repetative calls would be detrimental to
 *                  speed.
 *  dwrk double *   Double precision work array of length nu.
 *  nu   int        The number of unknowns and thus the dimensions of
 *                  lhs and rhs above.
 * Output:
 *  return int      0 - success, anything else means an error occurred.
 *                 -1 - singular matrix.
 */
int lu_solve(double **lhs, double *rhs, int *iwrk, double *dwrk, int nu)
{
/*
 * Decompose the input matrix into its LU representation.
 */
  if(lu_decomp(lhs, iwrk, dwrk, nu))
    return -1;
/*
 * Solve the equations via the LU decomposed matrix.
 */
  lu_backsub(lhs, rhs, iwrk, nu);
  return 0;
}

/*.......................................................................
 * Decompose a square 2D matrix into its Left-diagonalUpper-diagonal
 * version for use in lu_backsub(). LU decomposition is performed via
 * Crout's method along the lines described in Numerical Recipes in
 * FORTRAN, but note that unlike Numerical recipes, the input matrix
 * should be stored as an array of pointers to rows. (Numerical Recipes
 * in C doesn't take advantage of this possibility).
 *
 * Input/Output:
 *  lhs  double **  A SQUARE matrix organised as a C array of pointers
 *                  to row vectors, and dimensioned as lhs[nu][nu].
 *                  On output this contains the LU decomposition of
 *                  itself (row-swapped).
 * Output:
 *  indx    int *   The caller must supply an array of 'nu' elements.
 *                  On output this will record row-swaps.
 * Input:
 *  dwrk double *   Double precision work array of length nu.
 *  nu   int        Dimensions for 'lhs','indx' and 'dwrk'.
 * Output:
 *  return int      0 - success, anything else means an error occurred.
 *                 -1 - singular matrix.
 */
int lu_decomp(double **lhs, int *indx, double *dwrk, int nu)
{
  double absval;  /* The absolute value of a row element */
  double maxval;  /* Max absolute value */
  int maxpiv;     /* The row that holds the best pivot for the current row */
  int pivrow;     /* The row for which a pivot is sought */
  double dtmp;    /* Double precision temporary */
  double *ptmp;   /* Temporary pointer for swapping rows */
  double sum;
  int row;        /* Current row */
  int col;        /* Current column */
  int j;
/*
 * Determine the max (absolute) value on each input row. This will be used
 * for pivoting to determine the max normalised pivot on any row.
 * The scale factors will be stored in 'dwrk'.
 */
  for(row=0; row<nu; row++) {
    maxval = 0.0f;
    for(col=0; col<nu; col++) {
      absval = fabs(lhs[row][col]);
      if(absval > maxval) maxval=absval;
    };
    if(maxval == 0.0) {       /* Check for singular matrix */
      fprintf(stderr, "lu_decomp: Singular matrix.\n");
      return -1;
    };
    dwrk[row]=1.0/maxval;    /* Store the normalisation factor */
  };
/*
 * Loop over columns to implement Crout's method.
 */
  for(col=0; col<nu; col++) {
/*
 * Find the upper triangle coefficients (except for those on the diagonal)
 * for which sufficient lower triangle coefficients are known and store
 * them in the appropriate (no-longer required) element in 'lhs'.
 */
    for(row=0; row<col; row++) {
      sum=lhs[row][col];
      for(j=0; j<row; j++)
	sum -= lhs[row][j] * lhs[j][col];
      lhs[row][col]=sum;
    };
/*
 * Now find the lower triangle coefficients for which sufficient upper
 * triangle coefficients have been determined, and also the diagonal
 * for the upper triangle that couldn't be calculated above.
 */
    maxpiv=pivrow=col;
    maxval=0.0f;
    for(row=col; row<nu; row++) {
      sum=lhs[row][col];
      for(j=0; j<col; j++)
	sum -= lhs[row][j] * lhs[j][col];
      lhs[row][col]=sum;
/*
 * Now update search for a good pivot for row=column=col.
 */
      absval = dwrk[row]*fabs(sum);
      if(absval > maxval) {
	maxval=absval;
	maxpiv=row;
      };
    };
/*
 * Swap rows by swapping row pointers.
 */
    ptmp=lhs[pivrow]; lhs[pivrow] = lhs[maxpiv]; lhs[maxpiv] = ptmp;
/*
 * Swap the associated row normalization factors.
 */
    dtmp=dwrk[pivrow]; dwrk[pivrow]=dwrk[maxpiv]; dwrk[maxpiv]=dtmp;
/*
 * Record the row swap.
 */
    indx[pivrow]=maxpiv;
/*
 * Avoid singular matrix.
 */
    if(lhs[col][col]==0.0)
      lhs[col][col] = TINY;
/*
 * Divide the current column of the lower triangle by its pivot element.
 */
    dtmp = 1.0/lhs[col][col];
    for(j=col+1;j<nu;j++)
      lhs[j][col] *= dtmp;
  };
  return 0;
}

/*.......................................................................
 * Solve a set of linear simultaneous equations via LU decomposition
 * and backsubstitution, with the option of iterative improvement.
 *
 * Input:
 *  lhs  double **  A SQUARE matrix organised as a C array of pointers
 *                  to row vectors (dimensioned as lhs[nu][nu]).
 *                  This is the row-swapped, LU-decomposed left-hand-side
 *                  simultaneous equation matrix returned by lu_decomp().
 *  indx    int *   The row-swap index array returned by lu_decomp().
 *                  (dimensioned indx[nu]).
 *  nu   int        The dimension of the sent arrays.
 * Input/Output:
 *  rhs  double *   The original right-hand side vector of the simultaneous
 *                  equation matrix equation. (dimensioned rhs[nu]).
 *                  On output this contains the solution vector.
 */
void lu_backsub(double **lhs, double *rhs, int *indx, int nu)
{
  double sum;  /* Sum of determined parts of an equation */
  int row;     /* A row index */
  int col;     /* A column index */
  int piv;
  double dtmp;
/*
 * Row-swap 'rhs' to get a copy of 'rhs' in the same
 * order as the row-swapped 'lhs'.
 */
  for(row=0; row<nu; row++) {
    piv = indx[row];
    dtmp = rhs[row];
    rhs[row] = rhs[piv];
    rhs[piv] = dtmp;
  };
/*
 * Solve the upper triangle matrix equation via forward substitution.
 * Note that the diagonal elements of the upper triangle were not stored
 * in the LU decomposed array, but are known to be all 1.0 .
 */
  for(row=0;row<nu; row++) {
    sum=rhs[row];
    for(col=0; col<row; col++)
      sum -= lhs[row][col] * rhs[col];
    rhs[row] = sum;
  };
/*
 * Back-substitute for the final solutions.
 */
  for(row=nu-1;row>=0; row--) {
    sum=rhs[row];
    for(col=row+1; col<nu; col++)
      sum -= lhs[row][col] * rhs[col];
    rhs[row] = sum/lhs[row][row];
  };
  return;
}

/*.......................................................................
 * Take the LU decomposed version of a matrix and return the inverse of
 * the matrix.
 *
 * Input:
 *  lhs  double **  A SQUARE matrix organised as a C array of pointers
 *                  to row vectors (dimensioned as lhs[nu][nu]).
 *                  This is the row-swapped, LU-decomposed left-hand-side
 *                  simultaneous equation matrix returned by lu_decomp().
 *  indx    int *   The row-swap index array returned by lu_decomp().
 *                  (dimensioned indx[nu]).
 *  dwrk double *   Double precision work array of length nu.
 *  nu   int        The dimension of the sent arrays.
 * Input/Output:
 *  inv  double **  Send a square matrix, organised in the same manner as
 *                  'lhs'. This need not be initialized in any way.
 *                  On output it will contain the inverse of the matrix
 *                  of which 'lhs' is the LU decomposed version.
 *                  Note that 'inv' must not point at the same array as 'lhs'.
 * Output:
 *  return  int     0 - OK.
 *                  1 - Error.
 */
int lu_invert(double **lhs, double **inv, int *indx, double *dwrk, int nu)
{
  int row;  /* Row index */
  int col;  /* Column index */
/*
 * Clear the output matrix.
 */
  for(row=0; row<nu; row++) {
    for(col=0; col<nu; col++) {
      inv[row][col] = 0.0;
    };
  };
/*
 * Find the inverse of the matrix, one column at a time.
 */
  for(col=0; col<nu; col++) {
/*
 * Compose the associated column of the identity matrix in 'dwrk'.
 */
    for(row=0; row<nu; row++)
      dwrk[row] = 0.0;
    dwrk[col] = 1.0;
/*
 * Use the identity matrix column to solve for the associated inverse
 * matrix column.
 */
    lu_backsub(lhs, dwrk, indx, nu);
/*
 * Copy the inverse column into its position in the output inverse matrix.
 */
    for(row=0; row<nu; row++)
      inv[row][col] = dwrk[row];
  };
  return 0;
}
