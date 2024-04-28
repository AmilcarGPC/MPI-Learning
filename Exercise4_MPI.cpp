/*
PROGRAM: Calculate PI                                                            !
HOW TO RUN :                                            
$ mpicxx Exercise4_MPI.cpp -o run_Exercise4_MPI          
$ mpiexec -n 4 ./run_Exercise4_MPI
*/

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>
#include <time.h>

double f(double a) {
    return 4.0 / (1.0 + a * a);
}

int main(int argc, char* argv[]) {
    int n, myid, numprocs, i;
    double PI25DT = 3.141592653589793238462643;
    double mypi, pi, h, sum, x, cpu_time_used;
    clock_t start, end;
    double min_time, max_time, total_time;

    // MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // Processor 0 reads the number of intervals
    n = 1000000000;

    // MPI Broadcast n
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate interval size
    h = 1.0 / n;
    sum = 0.0;

    // Parallel calculation
    start = clock();
    for (i = myid + 1; i <= n; i += numprocs) {
        x = h * ((double)i - 0.5);
        sum += f(x);
    }

    mypi = h * sum;

    // MPI Reduction
    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    end = clock();

    // Calculate CPU time for each process
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    // MPI Reduction for CPU times
    MPI_Reduce(&cpu_time_used, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cpu_time_used, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cpu_time_used, &total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Processor 0 prints minimum, maximum, and average CPU times
    if (myid == 0) {
        printf("Min CPU time = %.4f sec\n", min_time);
        printf("Max CPU time = %.4f sec\n", max_time);
        printf("Average CPU time = %.4f sec\n", total_time / numprocs);
        printf("pi = %.15f, Error = %.15f\n", pi, fabs(pi - PI25DT));
    }

    // MPI Finalization
    MPI_Finalize();

    return 0;
}