/*
 File:   solve_flux_cuda.h
 Author: Kris Garrett
 Date:   April 19, 2013
*/

#ifndef __SOLVE_FLUX_CUDA_H
#define __SOLVE_FLUX_CUDA_H

extern "C"
void solveFluxInit_cuda(int nx, int ny, int nm, int nm2, int nq, double *w, double *p, 
                        double *xi, double *eta);

extern "C"
void solveFlux_cuda(int nx, int ny, int nm, int nq, double theta, double dx, double dy, 
                    double * alpha, double *flux);

#endif
