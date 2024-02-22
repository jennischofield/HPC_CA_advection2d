/*******************************************************************************
2D advection example program which advects a Gaussian u(x,y) at a fixed velocity



Outputs: initial.dat - initial values of u(x,y)
         final.dat   - final values of u(x,y)

         The output files have three columns: x, y, u

         Compile with: gcc -o advection2D -std=c99 advection2D.c -lm

Notes: The time step is calculated using the CFL condition

********************************************************************************/

/*********************************************************************
                     Include header files 
**********************************************************************/

#include <stdio.h>
#include <math.h>
#include <omp.h>


/*********************************************************************
                      Main function
**********************************************************************/

int main() {

    /* Grid properties */
    const int NX = 1000;    // Number of x points
    const int NY = 1000;    // Number of y points
    const float xmin = 0.0; // Minimum x value
    const float xmax = 30.0; // Maximum x value
    const float ymin = 0.0; // Minimum y value
    const float ymax = 30.0; // Maximum y value

    /* Parameters for the Gaussian initial conditions */
    const float x0 = 3.0;                    // Centre(x)
    const float y0 = 15.0;                    // Centre(y)
    const float sigmax = 1.0;               // Width(x)
    const float sigmay = 5.0;               // Width(y)
    const float sigmax2 = sigmax * sigmax; // Width(x) squared
    const float sigmay2 = sigmay * sigmay; // Width(y) squared

    /* Boundary conditions */
    const float bval_left = 0.0;    // Left boundary value
    const float bval_right = 0.0;   // Right boundary value
    const float bval_lower = 0.0;   // Lower boundary
    const float bval_upper = 0.0;   // Upper boundary

    /* Time stepping parameters */
    const float CFL = 0.9;   // CFL number
    const int nsteps = 800; // Number of time steps

    /* Velocity */
    const float velx = 1.0; // Velocity in x direction
    const float vely = 0.0; // Velocity in y direction

    /* Arrays to store variables. These have NX+2 elements
       to allow boundary values to be stored at both ends */
    float x[NX + 2];          // x-axis values
    float y[NX + 2];          // y-axis values
    float u[NX + 2][NY + 2];    // Array of u values
    float dudt[NX + 2][NY + 2]; // Rate of change of u


    float x2;   // x squared (used to calculate initial conditions)
    float y2;   // y squared (used to calculate initial conditions)

    /* Calculate distance between points */
    float dx = (xmax - xmin) / ((float) NX);
    float dy = (ymax - ymin) / ((float) NY);

    /* Calculate time step using the CFL condition */
    /* The fabs function gives the absolute value in case the velocity is -ve */
    float dt = CFL / ((fabs(velx) / dx) + (fabs(vely) / dy));

    /* 2.3 variables */
    double ustar = 0.2;
    double initialz = 1.0;
    double kappa = 0.41;
    float shear[NY + 2];
    /*** Report information about the calculation ***/
    printf("Grid spacing dx     = %g\n", dx);
    printf("Grid spacing dy     = %g\n", dy);
    printf("CFL number          = %g\n", CFL);
    printf("Time step           = %g\n", dt);
    printf("No. of time steps   = %d\n", nsteps);
    printf("End time            = %g\n", dt * (float) nsteps);
    printf("Distance advected x = %g\n", velx * dt * (float) nsteps);
    printf("Distance advected y = %g\n", vely * dt * (float) nsteps);

    /*** Place x points in the middle of the cell ***/
    /* LOOP 1 */
    // There are no dependencies within Loop 1, this loop can be
    // parallelised. Even though NX is a const, to compile correctly,
    // it must be declared as shared.
#pragma omp parallel for default (none) shared(x, dx)
    for (int i = 0; i < NX + 2; i++) {
        x[i] = ((float) i - 0.5) * dx;
    }

    /*** Place y points in the middle of the cell ***/
    /* LOOP 2 */
    /* There are no dependencies within Loop 2, this loop can be
     * parallelised. This loop has two shared variables - y and dy,
     * and no private variables. NY is a const and does not need to be
     * scoped according to the specification */
#pragma omp parallel for default (none) shared(y, dy)
    for (int j = 0; j < NY + 2; j++) {
        y[j] = ((float) j - 0.5) * dy;
    }

    /*Although not in the specification, this loop can also be parallelised.
     * it has no dependencies, and is just writing values to an array.
     * this loop determines the shear value and saves it to an array
     * for quick access in the time loop */
#pragma omp parallel for default (none) shared(y, initialz, shear, ustar, kappa)
    for (int k = 0; k < NY + 2; k++) {
        if (y[k] > initialz) {
            if (initialz != 0.0) {
                shear[k] = (double) (ustar / kappa) * log((y[k] / initialz));
            }
        } else if (y[k] <= initialz) {
            shear[k] = 0.0;
        }
    }
    /*Since we know that the initial list is sorted, we know that the value at NY + 1
    // is the greatest shear value. The shear value is then multiplied by velx
    // (in this case, velx is 1.0, but it gives more flexibility for experimentation
    // later on. We then determine the timestep using this maximum value.*/
    dt = CFL / ((fabs(shear[NY + 1] * velx) / dx) + (fabs(vely) / dy));

    /*** Set up Gaussian initial conditions ***/
    /* LOOP 3 */
    /*This loop can be parallelised, there are no data dependencies within this nested loop
     * By adding just the singular #pragma, it parallelises the outside loop, and each thread
     * evaluates the inner loop. It has the shared variables x,y, and u, and the private
     * variables x2 and y2, as well as constants NX, NY, x0, y0, sigmax2, and sigmay2 */
#pragma omp parallel for default (none) shared(x, y, u) private(x2, y2)
    for (int i = 0; i < NX + 2; i++) {
        for (int j = 0; j < NY + 2; j++) {
            x2 = (x[i] - x0) * (x[i] - x0);
            y2 = (y[j] - y0) * (y[j] - y0);
            u[i][j] = exp(-1.0 * ((x2 / (2.0 * sigmax2)) + (y2 / (2.0 * sigmay2))));
        }
    }

    /*** Write array of initial u values out to file ***/
    FILE *initialfile;
    initialfile = fopen("initial.dat", "w");
    /* LOOP 4 */
    /*This loop cannot be parallelised, it is an output dependency,
     * where each thread would have to write to initialfile.
     * However, if #pragma omp critical or #pragma omp ordered
     * was used, you could parallelise writing to a file, but that
     * isn't covered in the scope of this class */
    for (int i = 0; i < NX + 2; i++) {
        for (int j = 0; j < NY + 2; j++) {
            fprintf(initialfile, "%g %g %g\n", x[i], y[j], u[i][j]);
        }
    }
    fclose(initialfile);

    /*** Update solution by looping over time steps ***/
    /* LOOP 5 */
    /* No, this time loop can't be parallelised.
     * Internally, the array u is updated across multiple threads.
     * This is a flow dependency - each iteration of this loop requires the
     * previous loop to have completed and passed the updated value of u to it */
    for (int m = 0; m < nsteps; m++) {

        /*** Apply boundary conditions at u[0][:] and u[NX+1][:] ***/
        /* LOOP 6 */
        /*Yes, there are no dependencies. The value that changes from loop to loop is
         * independent of each thread. This loop contains shared variable u and
         * constants NY, bval_left, NX, and bval_right*/
#pragma omp parallel for default (none) shared(u)
        for (int j = 0; j < NY + 2; j++) {
            u[0][j] = bval_left;
            u[NX + 1][j] = bval_right;
        }

        /*** Apply boundary conditions at u[:][0] and u[:][NY+1] ***/
        /* LOOP 7 */
        /*Yes, there are no dependencies. The value that changes from loop to loop is
         * independent of each thread. This loop contains shared variable u and
         * constants NY, bval_lower, NX, and bval_upper*/
#pragma omp parallel for default (none)shared(u)
        for (int i = 0; i < NX + 2; i++) {
            u[i][0] = bval_lower;
            u[i][NY + 1] = bval_upper;
        }

        /*** Calculate rate of change of u using leftward difference ***/
        /* Loop over points in the domain but not boundary values */
        /* LOOP 8 */
        /*No, this is a flow dependency. To get u[i-1][j] or u[i][j-1]
         * you must calculate that before the current loop. */
        for (int i = 1; i < NX + 1; i++) {
            for (int j = 1; j < NY + 1; j++) {
                //Changed to multiply the velx variable by the vertical shear
                dudt[i][j] = (-(velx * shear[j]) * (u[i][j] - u[i - 1][j]) / dx
                              - vely * (u[i][j] - u[i][j - 1]) / dy);
            }
        }

        /*** Update u from t to t+dt ***/
        /* Loop over points in the domain but not boundary values */
        /* LOOP 9 */
        /* Yes, all values are internal to the loop. This loop has shared
         * variables u, dudt, and dt, and constants NX and NY*/
#pragma omp parallel for default (none) shared(u, dudt, dt)
        for (int i = 1; i < NX + 1; i++) {
            for (int j = 1; j < NY + 1; j++) {

                u[i][j] = u[i][j] + dudt[i][j] * dt;

            }
        }

    } // time loop

    /*** Write array of final u values out to file ***/

    FILE *finalfile;
    finalfile = fopen("final.dat", "w");
    /* LOOP 10 */
    /*This loop cannot be parallelised, it is an output dependency,
     * where each thread would have to write to finalfile.
     * However, if #pragma omp critical or #pragma omp ordered
     * was used, you could parallelise writing to a file, but that
     * isn't covered in the scope of this class */
    for (int i = 0; i < NX + 2; i++) {
        for (int j = 0; j < NY + 2; j++) {
            fprintf(finalfile, "%g %g %g\n", x[i], y[j], u[i][j]);
        }
    }
    fclose(finalfile);

    /*This loops writes out the averages for part 2.4 of the coursework to
     * a file for gnuplot to plot. */
    FILE *averagefile;
    averagefile = fopen("average.dat", "w");
    /*Even though not in the original template,
     * this loop cannot be parallelised, it is an output dependency,
     * where each thread would have to write to averagefile.
     * However, if #pragma omp critical or #pragma omp ordered
     * was used, you could parallelise writing to a file, but that
     * isn't covered in the scope of this class */
    for (int i = 0; i < NX + 2; i++) {
        double avgval = 0;
        for (int j = 0; j < NY + 2; j++) {
            avgval += u[i][j];
        }
        avgval /= NY + 1;
        fprintf(averagefile, "%g %g\n", x[i], avgval);
    }
    fclose(averagefile);
    return 0;
}

/* End of file ******************************************************/