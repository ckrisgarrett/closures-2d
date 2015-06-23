/*
 File:   fobj_cuda.cu
 Author: Kris Garrett
 Date:   June 7, 2012
*/

#ifndef __FOBJ_CUDA_H
#define __FOBJ_CUDA_H

extern "C"
int fobjUseCuda(int thread);

extern "C"
void fobjCudaInitSerial(int numCudaCards, int numThreadsPerCudaCard);

extern "C"
void fobjCudaInitParallel(int n, int nq, double *w, int thread);

extern "C"
double fobjMnCuda(int nm, int nq, double *alpha, double *u, double *p, int thread, 
                  double *g, double *h);

extern "C"
double fobjPPnCuda(int nm, int nq, double *alpha, double *u, double *w, double *p, 
                   double delta, int thread, double *g, double *h);

#endif
