/*
 File:   output.cpp
 Author: Kris Garrett
 Date:   Feb 01, 2013
*/

#include "kinetic_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdint.h>


/*
    Writes data to filename "out_<time>_<node>.sn".
    
    For each grid point, the vector of moments is written to the file.  
    The file consists of one long vector.
*/
void KineticSolver::outputData(double time)
{
    const int FILENAME_SIZE = 256;
    char filename[FILENAME_SIZE];
    FILE *file;
    int64_t temp;
    
    
    // Create and open file.
    snprintf(filename, FILENAME_SIZE, "out_%.3f_%d.sn", time, c_node);
    file = fopen(filename, "wb");
    
    // sizeX and sizeY of output grid
    temp = c_gX[2] - c_gX[1] + 1;
    fwrite(&temp, sizeof(int64_t), 1, file);
    temp = c_gY[2] - c_gY[1] + 1;
    fwrite(&temp, sizeof(int64_t), 1, file);
    
    // Domain bounds
    fwrite(&c_ax, sizeof(double), 1, file);
    fwrite(&c_ay, sizeof(double), 1, file);
    fwrite(&c_bx, sizeof(double), 1, file);
    fwrite(&c_by, sizeof(double), 1, file);
    
    // number of quadrature points and weights
    temp = c_numQuadPoints;
    fwrite(&temp, sizeof(int64_t), 1, file);
    fwrite(c_quadWeights, sizeof(double), c_numQuadPoints, file);
    
    // data
    for(int i = c_gX[1]; i <= c_gX[2]; i++)
    {
        for(int j = c_gY[1]; j <= c_gY[2]; j++)
        {
            fwrite(&c_kinetic[I3D(i,j,0)], sizeof(double), c_numQuadPoints, file);
        }
    }
    
    // Close file.
    fclose(file);
}

#ifdef USE_PAPI
void KineticSolver::papi_output() {
    printf("KineticSolver::update\n");
    papi_show_result(&c_update_info);
        printf("KineticSolver::solveFlux\n");
    papi_show_result(&c_flux_info);
}
#endif
