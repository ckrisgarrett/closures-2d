/*
 File:   utils.h
 Author: Kris Garrett
 Date:   July 26, 2012
*/

#ifndef __UTILS_H
#define __UTILS_H

#include <stdio.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

extern "C"
void utils_abort();

double utils_norm1(int n, double *A);
void utils_getGaussianWeightsAndNodes(int n, double *w, double *mu);
void utils_getLebedevWeightsAndNodes(int numPoints, double*, double*, double*);
int utils_numLebedevQuadPoints(int);

// Blas and lapack routine names.
extern "C"
void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x, 
            int *incX, double *beta, double *y, int *incY);
extern "C"
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
            double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc);
extern "C"
double dnrm2_(int *n, double *x, int *incX);
extern "C"
double dpotrf_(char *uplo, int *n, double *A, int *lda, int *info);
extern "C"
double dpocon_(char *uplo, int *n, double *A, int *lda, double *anorm, double *rcond, 
               double *work, int *iwork, int *info);
extern "C"
double dpotrs_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);
extern "C"
void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, 
            double *A, int *lda, double *B, int *ldb);
extern "C"
void dtrsv_(char *uplo, char *trans, char *diag, int *n, double *A, int *lda, double *x, int *incX);
extern "C"
void dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, 
            double *A, int *lda, double *B, int *ldb);
extern "C"
void dtrmv_(char *uplo, char *trans, char *diag, int *n, double *A, int *lda, double *x, int *incX);
extern "C"
void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, 
            double *x, int *incX, double *beta, double *y, int *incY);
extern "C"
double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
extern "C"
void dger_(int *m, int *n, double *alpha, double *x, int *incx, double *dy, int *incy, 
           double *A, int *lda);

#endif
