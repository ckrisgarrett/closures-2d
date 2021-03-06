/*
 File:   optcb.cpp
 Author: Kris Garrett
 Date:   February 26, 2012
*/

#include "opt.h"
#include "fobj.h"
#include "../../utils.h"
#include "linesearch.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_PAPI
#include "../../profiling.h"
#endif

static
void optScaled(int nm, int nq, double *u, double *P2M, double *alphaP, double *alphaM, 
               double *w, double *p, OPTIONS *options, OUTS *outs);



static void setIdentity(int n, double *A)
{
    memset(A, 0.0, n * n * sizeof(double));
    for(int i = 0; i < n; i++)
        A[i*n+i] = 1.0;
}

/*
    Calculate alpha for u = (1,0,0,...,0)^T
*/
static void isotropicAlpha(int n, double *alpha, int momentType, double delta)
{
    switch(momentType)
    {
        case MOMENT_TYPE_MN:
            alpha[0] = 2.0 * sqrt(M_PI) * log(0.5 / sqrt(M_PI));
            for(int i = 1; i < n; i++)
                alpha[i] = 0.0;
            break;
        case MOMENT_TYPE_PPN:
            alpha[0] = 1.0 - 4.0 * delta * M_PI;
            for(int i = 1; i < n; i++)
                alpha[i] = 0.0;
            break;
        default:
            printf("isotropicAlpha: moment type out of range.\n");
            utils_abort();
            break;
    }
}


/*
    Estimate an upper bound for the error.
*/
//static double modcfl(int n, double *da, int momentType)
//{
//    double da_1 = 0;
//    for(int i = 0; i < n; i++)
//    {
//        da_1 = da_1 + fabs(da[i]);
//    }
//    
//    int pnOrder = floor(sqrt(n)) + 1;
//    double m_inf = sqrt((2 * pnOrder + 1) / (2 * M_PI));
//    
//    switch(momentType)
//    {
//        case MOMENT_TYPE_MN:
//            return exp(2.0 * da_1 * m_inf);
//        case MOMENT_TYPE_PPN:
//            return 1.0 + 2.0 * da_1 * m_inf;
//        default:
//            printf("modcfl: moment type out of range.\n");
//            utils_abort();
//    }
//    
//    // Should never reach this point.
//    return 0.0;
//}


/*
    Scales u and calls the optimization routine.
*/
void optcb(int nm, int nq, double *u, double *P2M, double *alphaP, double *alphaM, 
           double *w, double *p, OPTIONS options, OUTS *outs)
{
    switch(options.momentType)
    {
        case MOMENT_TYPE_MN:
        {
            // Scale u and alpha.
            double u0 = u[0];
            double alphaScale = log(u0) * 2.0 * sqrt(M_PI);
            
            for(int i = 0; i < nm; i++)
                u[i] = u[i] / u0;
            
            if(u0 <= 0)
            {
                printf("u[0] < 0\n");
                utils_abort();
            }
            
            // Run optimization.
            optScaled(nm, nq, u, P2M, alphaP, alphaM, w, p, &options, outs);

            // Rescale u and alpha.
            for(int i = 0; i < nm; i++)
                u[i] = u[i] * u0;
            alphaM[0] = alphaM[0] + alphaScale;

            break;
        }
        case MOMENT_TYPE_PPN:
        {
            // Scale u and alpha.
            double u0 = u[0];
            options.deltaPPn = options.deltaPPn / (u0 * u0);
            
            for(int i = 0; i < nm; i++)
                u[i] = u[i] / u0;
            
            if(u0 <= 0)
            {
                printf("u[0] < 0\n");
                utils_abort();
            }
            
            // Run optimization.
            optScaled(nm, nq, u, P2M, alphaP, alphaM, w, p, &options, outs);
            
            // Rescale u and alpha.
            for(int i = 0; i < nm; i++)
            {
                alphaM[i] = alphaM[i] * u0;
                u[i] = u[i] * u0;
            }

            break;
        }
        default:
        {
            printf("Moment type not supported.\n");
            utils_abort();
        }
    }
}


/*
    Does the optimization for u = (1,...)^T
*/
void optScaled(int nm, int nq, double *u, double *P2M, double *alphaP, double *alphaM, 
               double *w, double *p, OPTIONS *options, OUTS *outs)
{
    #ifdef USE_PAPI
    profile_start_update("optcb.cpp: optScaled");
    #endif

    double f = 0.0;
    double *alphaP0 = (double*)malloc(nm * sizeof(double));
    double *dP = (double*)malloc(nm * sizeof(double));
    double *dM = (double*)malloc(nm * sizeof(double));
    double *cholesky = (double*)malloc(nm * nm * sizeof(double));
    double *gP = (double*)malloc(nm * sizeof(double));
    double *hP = (double*)malloc(nm * nm * sizeof(double));
    double *uP = (double*)malloc(nm * sizeof(double));
    double *gM = (double*)malloc(nm * sizeof(double));
    double *P2Min = (double*)malloc(nm * nm * sizeof(double));
    double *pP = (double*)malloc(nq * nm * sizeof(double));
    int inc1 = 1;
    char lChar = 'L';
    char rChar = 'R';
    char tChar = 'T';
    char nChar = 'N';
    double one = 1.0;
    
    memcpy(alphaP0, alphaP, nm * sizeof(double));
    memcpy(P2Min, P2M, nm * nm * sizeof(double));
    
    // pP holds p in the P basis
    memcpy(pP, p, nq * nm * sizeof(double));
    dtrsm_(&lChar, &lChar, &nChar, &nChar, &nm, &nq, &one, P2Min, &nm, pP, &nm);
    
    int maxIter = options->maxIter;
    double tolRel = options->tolRel;
    //double tolGamma = options->tolGamma;
    double r = 0.0;
    double condHMax = options->condHMax;
    int regIndex = 0;
    
    // uP <- P2M \ u
    memcpy(uP, u, nm * sizeof(double));
    dtrsv_(&lChar, &nChar, &nChar, &nm, P2M, &nm, uP, &inc1);
    
    
    // tolerance relative to |u|
    double rtol = dnrm2_(&nm, u, &inc1) * tolRel;
    
    
    // initialize the iteration counter
    int iter = 0;
    int totalIter = 0; // Counts iterations including regularization.
    
    // initialize a variable to hold the number of steps we took get get gamma
    // below its tolerance after the norm of the gradient satisfies its
    // tolerance.  (typically, the norm of the gradient will fall below its
    // tolerance much more quickly than gamma will.)
    //int iterGamma = 0;

    // initialize the backtracking factor
    double t = 1;
    //double gamma, err;
    double err = 0.0;
    //gamma = err = 0;      // To stop compiler warnings.
    
    
    // To see if we hit max regularization once already.
    bool hitMaxRegularization = false;
    
    
    // Iterate to find alpha.
    while(true)
    {
        // increment the iteration counters
        iter = iter + 1;
        totalIter = totalIter + 1;
        
        // Find f, g, and h.
        switch(options->momentType)
        {
            case MOMENT_TYPE_MN:
                f = fobjMn(nm, nm, nq, alphaP, uP, w, pP, gP, hP, false);
                break;
            case MOMENT_TYPE_PPN:
                f = fobjPPn(nm, nm, nq, alphaP, uP, w, pP, options->deltaPPn, gP, hP, false);
        }

        
        // Check quadrature and find Newton direction with Cholesky factorization.
        // ALSO NEED TO CHECK CONDITION NUMBER OF H.
        memcpy(cholesky, hP, nm * nm * sizeof(double));
        int error;
        dpotrf_(&lChar, &nm, cholesky, &nm, &error);
        double norm = utils_norm1(nm, hP);
        double condest = 0.0;
        int error2;
        double *work = (double*)malloc(3*nm*sizeof(double));
        int *iwork = (int*)malloc(nm*sizeof(int));
        dpocon_(&lChar, &nm, cholesky, &nm, &norm, &condest, work, iwork, &error2);
        free(work);
        free(iwork);
        condest = 1.0 / condest;
        
        
        // Initialize to isotropic alpha and no change of basis if 
        // H is too bad on first iteration.
        if(totalIter == 1 && (error != 0 || condest > condHMax))
        {
            // Reset Variables with standard basis.
            isotropicAlpha(nm, alphaP0, options->momentType, options->deltaPPn);
            setIdentity(nm, P2M);
            setIdentity(nm, P2Min);
            
            memcpy(alphaP, alphaP0, nm * sizeof(double));
            memcpy(uP, u, nm * sizeof(double));
            memcpy(pP, p, nq * nm * sizeof(double));
            
            rtol = dnrm2_(&nm, u, &inc1) * tolRel;
            iter = 0;
            //iterGamma = 0;
            t = 1;
            
            continue;
        }

        
        // Regularize criteria.
        if(iter >= maxIter || error != 0 || t == 0 || condest > condHMax)
        {
            // REGULARIZE ...
            if(regIndex < NUM_REGULARIZATIONS - 1)
            {
                regIndex++;
                r = REGULARIZATIONS[regIndex];
                
                uP[0] = 1.0;
                for(int i = 1; i < nm; i++)
                    uP[i] = u[i] * (1.0 - r);
            }
            else
            {
                if(hitMaxRegularization)
                {
                    // This should never happen.
                    printf("Max regularization achieved.\n");
                    utils_abort();
                }
                isotropicAlpha(nm, alphaP0, options->momentType, options->deltaPPn);
                setIdentity(nm, P2Min);
                memcpy(uP, u, nm * sizeof(double));
                regIndex = 0;
                hitMaxRegularization = true;
                maxIter = options->maxIter;
            }
            
            
            rtol = dnrm2_(&nm, uP, &inc1) * tolRel;
            memcpy(alphaP, alphaP0, nm * sizeof(double));
            memcpy(P2M, P2Min, nm * nm * sizeof(double));
            dtrsv_(&lChar, &nChar, &nChar, &nm, P2M, &nm, uP, &inc1);
            
            memcpy(pP, p, nq * nm * sizeof(double));
            dtrsm_(&lChar, &lChar, &nChar, &nChar, &nm, &nq, &one, P2Min, &nm, pP, &nm);
            
            iter = 0;
            //iterGamma = 0;
            t = 1;
            if(regIndex == NUM_REGULARIZATIONS - 1)
                maxIter = 10 * maxIter;
            
            continue;
        }
        
        
        // H = R^T * R, where R is the upper triangular part of cholesky
        dtrsv_(&lChar, &nChar, &nChar, &nm, cholesky, &nm, gP, &inc1);
        // change the basis of the vector holding the basis polynomials
        dtrsm_(&lChar, &lChar, &nChar, &nChar, &nm, &nq, &one, cholesky, &nm, pP, &nm);
        // update P2M <- P2M * R
        dtrmm_(&rChar, &lChar, &nChar, &nChar, &nm, &nm, &one, cholesky, &nm, P2M, &nm);
        // change alpha's basis: alphaP <- R * alphaP
        dtrmv_(&lChar, &tChar, &nChar, &nm, cholesky, &nm, alphaP, &inc1);
        // change the basis of v(r): v_P(r) <- R^T \ v_P(r)
        dtrsv_(&lChar, &nChar, &nChar, &nm, cholesky, &nm, uP, &inc1);
        
        
        // Solve for the newton direction d = -g in the new basis
        for(int i = 0; i < nm; i++)
            dP[i] = -gP[i];


        // the error as measured by the norm of the gradient in the legendre basis
        memcpy(gM, gP, nm * sizeof(double));
        dtrmv_(&lChar, &nChar, &nChar, &nm, P2M, &nm, gM, &inc1);
        err = dnrm2_(&nm, gM, &inc1);
        
        
        // now check the stopping criteria
        if(err <= rtol)
        {
            // gamma, the supplementary stopping criterion
            memcpy(dM, dP, nm * sizeof(double));
            dtrsv_(&lChar, &tChar, &nChar, &nm, P2M, &nm, dM, &inc1);
            //gamma = modcfl(nm, dM, options->momentType);
            //if(gamma <= tolGamma)
            {
                // problem solved!
                memcpy(alphaM, alphaP, nm * sizeof(double));
                dtrsv_(&lChar, &tChar, &nChar, &nm, P2M, &nm, alphaM, &inc1);
                
                u[0] = 1.0;
                for(int i = 1; i < nm; i++)
                    u[i] = u[i] * (1.0 - r);
                
                break;
            }
        }
//        else if(err <= rtol && gamma > tolGamma)
//        {
//            // this is to keep track of how many extra iterations we took just for gamma
//            iterGamma++;
//        }
        
        
        // if you get to this point, that means you didn't satisfy the stopping
        // criteria and you didn't need to regularize more than you already are
        // so, it's time to take a step
        t = linesearch(nm, nm, nq, alphaP, dP, f, gP, uP, w, pP, 
            options->deltaPPn, options->momentType);
    }
    
    
    if(outs != NULL)
    {
        outs->iter = totalIter;
//        outs->iterGamma = iterGamma;
//        if(outs->iterGamma < 0)
//            outs->iterGamma = 0;
//        outs->gamma = gamma;
        outs->normG = err;
        outs->r = r;
    }
    
    #ifdef USE_PAPI
    profile_finish_update("optcb.cpp: optScaled");
    #endif
    
    free(pP);
    free(dM);
    free(dP);
    free(cholesky);
    free(gP);
    free(gM);
    free(hP);
    free(uP);
    free(alphaP0);
    free(P2Min);
}

