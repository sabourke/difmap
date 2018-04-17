#ifndef matinv_h
#define matinv_h

int gj_invert(double **lhs, int *work1, int *work2, int *work3, int nu);
int gj_solve(double **lhs, double *rhs, int *work1, double **work2, int nu);

int lu_solve(double **lhs, double *rhs, int *iwrk, double *dwrk, int nu);
int lu_decomp(double **lhs, int *indx, double *dwrk, int nu);
void lu_backsub(double **lhs, double *rhs, int *indx, int nu);
int lu_invert(double **lhs, double **inv, int *indx, double *dwrk, int nu);

#endif
