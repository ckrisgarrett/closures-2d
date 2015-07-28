/*
 File:   fobj.cpp
 Author: Kris Garrett
 Date:   May 25, 2012
*/

#include "fobj.h"
#include <math.h>
#include <stdlib.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif


HStruct *g_hStruct;     // Global variable for Gaunt coefficients.
static double *gExtendedForAllThreads;


/*
    Initializes static variables.
*/
void initFobj(int numOmpThreads, int numManyMoments)
{
    gExtendedForAllThreads = new double[numOmpThreads * numManyMoments];
}


/*
    Solves for the Hessian for the Mn algorithm.
*/
//__attribute__ ((noinline))  // This is here to aid in profiling.
static
void solveHMn(int nm, int np2, int nq, double *alpha, double *w, double *p, 
              double *h, bool useGaunt)
{
    int threadId = 0;
    double *gExtended;
    
    
    // Get OpenMP thread index.
    #ifdef USE_OPENMP
    threadId = omp_get_thread_num();
    #endif
    
    
    // Zero variables.
    gExtended = &gExtendedForAllThreads[threadId*np2];
    for(int k = 0; k < np2; k++)
        gExtended[k] = 0.0;
    
    for(int i = 0; i < nm*nm; i++)
        h[i] = 0.0;
    
    
    // Integrate.
    for(int q = 0; q < nq; q++)
    {
        // Compute temp = exp(alpha^T * p).
        double temp = 0;
        for(int i = 0; i < nm; i++)
        {
            temp = temp + alpha[i] * p[q*np2+i];
        }
        temp = exp(temp);

        
        // Compute H.
        if(useGaunt)
        {
            for(int k = 0; k < np2; k++)
            {
                gExtended[k] += w[q] * p[q*np2+k] * temp;
            }
        }
        else
        {
            for(int i = 0; i < nm; i++)
            {
                double pi = p[q*np2+i];
                for(int j = 0; j <= i; j++)
                {
                    double pj = p[q*np2+j];
                    h[i*nm+j] += w[q] * pi * pj * temp;
                }
            }
        }
    }
    
    
    // Finish the computations.
    if(useGaunt)
    {
        for(int i = 0; i < nm; i++)
        {
            for(int j = i; j < nm; j++)
            {
                double h_ij = 0.0;
                int index = i * nm + j;
                for(int k = 0; k < g_hStruct[index].numMoments; k++)
                {
                    h_ij += g_hStruct[index].s[k].value * gExtended[g_hStruct[index].s[k].moment];
                }
                h[i*nm+j] = h[j*nm+i] = h_ij;
            }
        }
    }
    else
    {
        for(int i = 0; i < nm; i++)
        {
            for(int j = i+1; j < nm; j++)
            {
                h[i*nm+j] = h[j*nm+i];
            }
        }
    }
}


/*
    Solves for the Hessian for the Mn algorithm.
*/
//__attribute__ ((noinline))  // This is here to aid in profiling.
static
void solveHPPn(int nm, int np2, int nq, double *alpha, double *w, double *p, double delta, 
               double *h, bool useGaunt)
{
    int threadId = 0;
    double *gExtended;
    
    
    // Get OpenMP thread index.
    #ifdef USE_OPENMP
    threadId = omp_get_thread_num();
    #endif
    
    
    // Zero variables.
    gExtended = &gExtendedForAllThreads[threadId*np2];
    for(int k = 0; k < np2; k++)
        gExtended[k] = 0.0;
    
    for(int i = 0; i < nm*nm; i++)
        h[i] = 0.0;
    
    
    // Integrate.
    for(int q = 0; q < nq; q++)
    {
        // Compute temp = 0.5 * alpha^T * p.
        double temp = 0;
        for(int i = 0; i < nm; i++)
        {
            temp = temp + alpha[i] * p[q*np2+i];
        }
        temp = 0.5 * temp;
        
        // Compute temp2 = sqrt((alpha^t * p)^2 + delta).
        double temp2 = sqrt(temp * temp + delta);
        double temp3 = temp + temp2;
        if(temp < 0)
            temp3 = -delta / (temp - sqrt(temp*temp + delta));
        double temp4 = 0.5 * w[q] * (temp3 / temp2);

        
        // Compute H.
        if(useGaunt)
        {
            for(int k = 0; k < np2; k++)
            {
                gExtended[k] += p[q*np2+k] * temp4;
            }
        }
        else
        {
            for(int i = 0; i < nm; i++)
            {
                double pi = p[q*np2+i];
                for(int j = 0; j <= i; j++)
                {
                    double pj = p[q*np2+j];
                    h[i*nm+j] += pi * pj * temp4;
                }
            }
        }
    }
    
    
    // Finish the computations.
    if(useGaunt)
    {
        for(int i = 0; i < nm; i++)
        {
            for(int j = i; j < nm; j++)
            {
                double h_ij = 0.0;
                int index = i * nm + j;
                for(int k = 0; k < g_hStruct[index].numMoments; k++)
                {
                    h_ij += g_hStruct[index].s[k].value * gExtended[g_hStruct[index].s[k].moment];
                }
                h[i*nm+j] = h[j*nm+i] = h_ij;
            }
        }
    }
    else
    {
        for(int i = 0; i < nm; i++)
        {
            for(int j = i+1; j < nm; j++)
            {
                h[i*nm+j] = h[j*nm+i];
            }
        }
    }
}


/*
    Evaluates f = integral(e^(alpha^T p)) - alpha^T u
              g = integral(p e^(alpha^T p)) - u
              h = integral(p p^T e^(alpha^T p))
    
    Variables:
    alpha  Vector (size number of moments)
    u      Vector (size number of moments)
    w      Vector (size number of quadrature points)
    p      Matrix of moments evaluated at quadrature points (size quadrature by moments)
    g      Gradient of f with respect to alpha (size moments)
    h      Hessian of f with respect to alpha (size moments by moments)

    Notes:
    If h != NULL, then g is assumed not NULL.
*/
double fobjMn(int nm, int np2, int nq, double *alpha, double *u, double *w, double *p, 
              double *g, double *h, bool useGaunt)
{
    double f = 0;
    
    
    // If H != NULL, then get all the information for f and g from H.
    if(h != NULL)
    {
        solveHMn(nm, np2, nq, alpha, w, p, h, useGaunt);
        
        if(g != NULL)
        {
            for(int i = 0; i < nm; i++)
            {
                g[i] = h[i] / p[0] - u[i];
            }
        }
        
        f = h[0] / (p[0] * p[0]);
        for(int i = 0; i < nm; i++)
            f = f - alpha[i] * u[i];
        
        return f;
    }
    
    
    // If H == NULL, then find g and f.
    if(g != NULL)
    {
        for(int i = 0; i < nm; i++)
            g[i] = 0.0;
    }
    
    for(int q = 0; q < nq; q++)
    {
        // Compute temp = exp(alpha^T * p).
        double temp = 0;
        for(int i = 0; i < nm; i++)
        {
            temp = temp + alpha[i] * p[q*np2+i];
        }
        temp = exp(temp);
        
        
        if(g != NULL)
        {
            for(int i = 0; i < nm; i++)
            {
                g[i] += w[q] * p[q*np2+i] * temp;
            }
        }

        f = f + w[q] * temp;
    }
    
    
    if(g != NULL)
    {
        for(int i = 0; i < nm; i++)
            g[i] = g[i] - u[i];
    }
    for(int i = 0; i < nm; i++)
    {
        f = f - alpha[i] * u[i];
    }

    return f;
}


/*
    Evaluates f = integral(0.5 ap (ap/2 + 0.5*sqrt(ap^2 + 4delta) + 0.5 delta) - delta + 
                           delta log(ap/2 + 0.5 sqrt(ap^2 + 4 delta))) - alpha^T u
              (ap = alpha^T p for the 'f' formula.)
    
    Variables:
    alpha  Vector (size number of moments)
    u      Vector (size number of moments)
    w      Vector (size number of quadrature points)
    p      Matrix of moments evaluated at quadrature points (size quadrature by moments)
    delta  The delta factor for the PPn method
    g      Gradient of f with respect to alpha (size moments)
    h      Hessian of f with respect to alpha (size moments by moments)

    Notes:
    If h != NULL, then g is assumed not NULL.
*/
double fobjPPn(int nm, int np2, int nq, double *alpha, double *u, double *w, double *p, double delta, 
               double *g, double *h, bool useGaunt)
{
    double f = 0;
    if(g != NULL)
    {
        for(int i = 0; i < nm; i++)
            g[i] = 0.0;
    }
    if(h != NULL)
    {
        for(int i = 0; i < nm*nm; i++)
            h[i] = 0.0;
    }
    

    // Integrate.
    for(int q = 0; q < nq; q++)
    {
        // Compute temp = 0.5 * alpha^T * p.
        double temp = 0;
        for(int i = 0; i < nm; i++)
        {
            temp = temp + alpha[i] * p[q*np2+i];
        }
        temp = 0.5 * temp;
        
        // Compute temp2 = sqrt((alpha^t * p)^2 + delta).
        double temp2 = sqrt(temp * temp + delta);
        double temp3 = temp + temp2;
        if(temp < 0)
            temp3 = -delta / (temp - sqrt(temp * temp + delta));
        
        // Compute f.
        f = f + w[q] * (temp * (temp3) - 0.5 * delta + delta * log(temp3));

        // Compute g.
        if(g != NULL)
        {
            double temp4 = w[q] * (temp3);
            for(int i = 0; i < nm; i++)
            {
                g[i] += p[q*np2+i] * temp4;
            }
        }
    }
    
    
    // Compute H.
    if(h != NULL)
    {
        solveHPPn(nm, np2, nq, alpha, w, p, delta, h, useGaunt);
    }


    // Finish the computations not requiring integration.
    for(int i = 0; i < nm; i++)
    {
        f = f - alpha[i] * u[i];

        if(g != NULL)
            g[i] = g[i] - u[i];
    }

    return f;
}



