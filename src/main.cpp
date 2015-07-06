/*
 File:   main.cpp
 Author: Kris Garrett
 Date:   July 26, 2012
*/

#include <stdio.h>
#include "utils.h"
#include "input_deck_reader.h"
#include "moment/moment_solver.h"
#include "kinetic/kinetic_solver.h"
#include "momopt/momopt_solver.h"
#include "dn/dn_solver.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

/*
    Checks the status of input parameters.
*/
static void checkInput(bool isOk, int lineNumber)
{
    if(!isOk) {
        printf("main.cpp:input error at line %d\n", lineNumber);
        utils_abort();
    }
}

/*
    Start the program.
*/
int main(int argc, char **argv)
{
    InputDeckReader inputDeckReader;
    char solverType[100];
    int numCellsX, numCellsY, numNodesX, numNodesY, initCond;
    double ax, ay, bx, by, tFinal, outDeltaT, gaussianSigma, initfloor, sigma;

    // Initialize MPI.
    #ifdef USE_MPI
    MPI_Init(&argc, &argv);
    #endif
    

    // Get the node index and number of nodes.
    int node = 0;
    int numNodes = 1;
    #ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
    #endif
    

    // Read input.deck
    inputDeckReader.readInputDeck("input.deck");
    checkInput(inputDeckReader.getValue("SOLVER", solverType), __LINE__);
    checkInput(inputDeckReader.getValue("NUM_CELLS_X", &numCellsX), __LINE__);
    checkInput(inputDeckReader.getValue("NUM_CELLS_Y", &numCellsY), __LINE__);
    checkInput(inputDeckReader.getValue("NUM_MPI_PARTITIONS_X", &numNodesX), __LINE__);
    checkInput(inputDeckReader.getValue("NUM_MPI_PARTITIONS_Y", &numNodesY), __LINE__);
    checkInput(inputDeckReader.getValue("A_X", &ax), __LINE__);
    checkInput(inputDeckReader.getValue("B_X", &bx), __LINE__);
    checkInput(inputDeckReader.getValue("A_Y", &ay), __LINE__);
    checkInput(inputDeckReader.getValue("B_Y", &by), __LINE__);
    checkInput(inputDeckReader.getValue("T_FINAL", &tFinal), __LINE__);
    checkInput(inputDeckReader.getValue("OUT_DELTA_T", &outDeltaT), __LINE__);
    checkInput(inputDeckReader.getValue("GAUSSIAN_SIGMA", &gaussianSigma), __LINE__);
    checkInput(inputDeckReader.getValue("FLOOR", &initfloor), __LINE__);
    checkInput(inputDeckReader.getValue("INIT_COND", &initCond), __LINE__);
    checkInput(inputDeckReader.getValue("SIGMA", &sigma), __LINE__);


    // Print out input.deck
    if(node == 0)
    {
        inputDeckReader.print();
    }


    // Create the solver.
    Solver *solver = 0;
    if(strcmp(solverType, "kinetic") == 0)
        solver = new KineticSolver();
    else if(strcmp(solverType, "moment") == 0)
        solver = new MomentSolver();
    else if(strcmp(solverType, "momopt") == 0)
        solver = new MomOptSolver();
    else if(strcmp(solverType, "dn") == 0)
        solver = new DnSolver();
    else {
        printf("Solver type does not exist\n.");
        utils_abort();
    }
    solver->c_gaussianSigma = gaussianSigma;
    solver->c_floor = initfloor;
    solver->c_initCond = initCond;
    solver->c_inputDeckReader = inputDeckReader;
    
    // Set the number of openmp threads.
    // This must be after readInputDeck and before init.
    #ifdef USE_OPENMP
    #pragma omp parallel
    {
        if(omp_get_thread_num() == 0)
        {
            solver->c_numOmpThreads = omp_get_num_threads();
            omp_set_num_threads(solver->c_numOmpThreads);
        }
    }
    #else
    solver->c_numOmpThreads = 1;
    #endif
    
    if(node == 0)
        printf("OpenMP Threads: %d\n", solver->c_numOmpThreads);


    // Set the size of each node's domain in the x/y coordinates.
    double dx = (bx-ax) / numCellsX;
    double dy = (by-ay) / numCellsY;
    solver->c_sizeX = numCellsX;
    solver->c_sizeY = numCellsY;
    solver->c_node = node;
    solver->c_numNodes = numNodes;
    solver->c_nodeX = 0;
    solver->c_nodeY = 0;
    solver->c_numNodesX = 1;
    solver->c_numNodesY = 1;
    solver->c_ax = ax;
    solver->c_ay = ay;
    solver->c_bx = bx;
    solver->c_by = by;
    solver->c_initX = ax;
    solver->c_initY = ay;

    #ifdef USE_MPI
    solver->c_numNodesX = numNodesX;
    solver->c_numNodesY = numNodesY;
    if(numCellsX < numNodesX || numCellsY < numNodesY)
    {
        printf("Init: Number of partitions larger than number of cells.\n");
        utils_abort();
    }
    if(numNodesX * numNodesY != numNodes)
    {
        printf("Init: Number of Nodes != numPartitionsX * numPartitionsY.\n");
        utils_abort();
    }
    solver->c_nodeX = node / numNodesY;
    solver->c_nodeY = node % numNodesY;

    solver->c_sizeX = solver->c_sizeX / numNodesX;
    if(solver->c_nodeX < numCellsX % numNodesX)
    {
        solver->c_sizeX++;
        solver->c_initX = ax + dx * solver->c_sizeX * solver->c_nodeX;
    }
    else
    {
        solver->c_initX = bx - dx * solver->c_sizeX * (solver->c_numNodesX - solver->c_nodeX);
    }
    solver->c_sizeY = solver->c_sizeY / numNodesY;
    if(solver->c_nodeY < numCellsY % numNodesY)
    {
        solver->c_sizeY++;
        solver->c_initY = ay + dy * solver->c_sizeY * solver->c_nodeY;
    }
    else
    {
        solver->c_initY = by - dy * solver->c_sizeY * (solver->c_numNodesY - solver->c_nodeY);
    }
    #endif
    

    // Initialize data.
    solver->initializeGrid(solver->getNumGhostCells(), dx, dy, sigma);
    double max_dt = solver->init(dx, dy);
    if(node == 0)
        printf("Init done.\n");
    
    
    // Start loop in time.
    double t = 0;
    double dt = max_dt;
    solver->outputData(t);
    for(double tOut = MIN(t + outDeltaT, tFinal); tOut <= tFinal; 
        tOut = MIN(tOut + outDeltaT, tFinal))
    {
        dt = max_dt;
        while(t < tOut)
        {
            t += dt;
            if(t > tOut)
            {
                dt = tOut - (t - dt);
                t = tOut;
            }
            
            solver->update(dt, dx, dy);
            
            if(node == 0 && (t*100 - floor(t*100)) / 100.0 < dt)
                printf("dt = %f   t = %f\n", dt, t);
        }
        
        solver->outputData(t);
        if(tOut > tFinal - 1e-12 * tFinal)
            break;
    }

    
    // Print time taken by the program.
    #ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    // End MPI.
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    
    return 0;
}


/*
    Initializes the grid.
*/
void Solver::initializeGrid(int numGhostCells, double dx, double dy, double sigma)
{
    // X and Y indices for where ghost cells start, real cells start,
    // real cells end, and ghost cells end.
    c_gX[0] = 0;
    c_gX[1] = numGhostCells;
    c_gX[2] = numGhostCells + c_sizeX - 1;
    c_gX[3] = 2 * numGhostCells + c_sizeX - 1;
    
    c_gY[0] = 0;
    c_gY[1] = numGhostCells;
    c_gY[2] = numGhostCells + c_sizeY - 1;
    c_gY[3] = 2 * numGhostCells + c_sizeY - 1;
    
    int totalSizeX = c_gX[3] - c_gX[0] + 1;
    int totalSizeY = c_gY[3] - c_gY[0] + 1;
    
    
    // Allocate initial grid.
    c_initialGrid = (double*)malloc(totalSizeX * totalSizeY * sizeof(double));
    c_sigmaT      = (double*)malloc(totalSizeX * totalSizeY * sizeof(double));
    c_sigmaS      = (double*)malloc(totalSizeX * totalSizeY * sizeof(double));
    
    
    // Initial conditions
    for(int i = c_gX[0]; i <= c_gX[3]; i++)
    {
        for(int j = c_gY[0]; j <= c_gY[3]; j++)
        {
            double x_i = c_initX + (i-c_gX[1]) * dx + dx/2.0;
            double y_j = c_initY + (j-c_gY[1]) * dy + dy/2.0;
                        
            
            // Boundary
            if(i < c_gX[1] || i > c_gX[2] || 
               j < c_gY[1] || j > c_gY[2])
            {
                c_initialGrid[I2D(i,j)] = c_floor;
                c_sigmaT[I2D(i,j)] = sigma;
                c_sigmaS[I2D(i,j)] = sigma;
            }
            
            
            // Delta Initial Condition.
            else if(c_initCond == INITCOND_LINESOURCE && c_gaussianSigma == 0.0)
            {
                c_initialGrid[I2D(i,j)] = c_floor;
                if(c_nodeX == (c_numNodesX - 1) / 2 &&
                   c_nodeY == (c_numNodesY - 1) / 2)
                {
                    if((c_numNodesX % 2 != 0 && i == (c_sizeX + numGhostCells) / 2) ||
                       (c_numNodesX % 2 == 0 && i == c_gX[2]))
                    {
                        if((c_numNodesY % 2 != 0 && j == (c_sizeY + numGhostCells) / 2) ||
                           (c_numNodesY % 2 == 0 && j == c_gY[2]))
                        {
                            c_initialGrid[I2D(i,j)] = 1.0 / (dx * dy);
                        }
                    }
                }
                
                c_sigmaT[I2D(i,j)] = sigma;
                c_sigmaS[I2D(i,j)] = sigma;
            }
            
            
            // Gaussian Initial Condition
            else if(c_initCond == INITCOND_LINESOURCE)
            {
                c_initialGrid[I2D(i,j)] = c_floor;
                double s2 = c_gaussianSigma * c_gaussianSigma;
                double gaussianFactor = 1.0 / (2.0 * M_PI * s2) * 
                    exp(-(x_i * x_i + y_j * y_j) / (2.0 * s2));
                c_initialGrid[I2D(i,j)] = MAX(gaussianFactor, c_floor);
                
                c_sigmaT[I2D(i,j)] = sigma;
                c_sigmaS[I2D(i,j)] = sigma;
            }
            
            
            // Lattice Initial Condition
            else if(c_initCond == INITCOND_LATTICE)
            {
                c_initialGrid[I2D(i,j)] = c_floor;
                c_sigmaT[I2D(i,j)] = 1.0;
                c_sigmaS[I2D(i,j)] = 1.0;
                
                if( (fabs(x_i - 2.0) < 0.5 && fabs(y_j - 2.0) < 0.5) || 
                    (fabs(x_i - 2.0) < 0.5 && fabs(y_j - 0.0) < 0.5) || 
                    (fabs(x_i - 2.0) < 0.5 && fabs(y_j + 2.0) < 0.5) || 
                    (fabs(x_i - 1.0) < 0.5 && fabs(y_j - 1.0) < 0.5) || 
                    (fabs(x_i - 1.0) < 0.5 && fabs(y_j + 1.0) < 0.5) || 
                    (fabs(x_i - 0.0) < 0.5 && fabs(y_j + 2.0) < 0.5) || 
                    (fabs(x_i + 1.0) < 0.5 && fabs(y_j - 1.0) < 0.5) || 
                    (fabs(x_i + 1.0) < 0.5 && fabs(y_j + 1.0) < 0.5) || 
                    (fabs(x_i + 2.0) < 0.5 && fabs(y_j - 2.0) < 0.5) || 
                    (fabs(x_i + 2.0) < 0.5 && fabs(y_j - 0.0) < 0.5) || 
                    (fabs(x_i + 2.0) < 0.5 && fabs(y_j + 2.0) < 0.5) )
                {
                    c_sigmaT[I2D(i,j)] = 10.0;
                    c_sigmaS[I2D(i,j)] = 0.0;
                }
            }
            
            
            // Smooth Periodic Condition
            else if(c_initCond == INITCOND_PERIODIC)
            {
                c_sigmaT[I2D(i,j)] = sigma;
                c_sigmaS[I2D(i,j)] = sigma;
                
                c_initialGrid[I2D(i,j)] = 1 + sin(2 * M_PI * x_i) * cos(2 * M_PI * y_j);
            }
        }
    }
}
