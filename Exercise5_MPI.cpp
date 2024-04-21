/*
PROGRAM: Calculate integral                                                          !
HOW TO RUN :                                            
$ mpicxx Exercise5_MPI.cpp -o Exercise5_MPI          
$ mpiexec -n 4 ./Exercise5_MPI
*/

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

double f(double x, double y) {
    return 1 + 6 * x * y * y;
}

int main(int argc, char* argv[]) {
    int n, m, myid, numprocs, i, j;
    double result, dx, dy, sum,gsum, x, y, cpu_time_used;
    clock_t start, end;
    double min_time, max_time, total_time;

    // MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // Processor 0 reads the number of intervals
    n = 10000;
    m = 10000;
    double a = 0.0, b = 2.0, c = -1.0, d = 1.0;

    // MPI Broadcast n
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate interval size
    dx = (b - a) / m;
    dy = (d - c) / n;
    sum = 0.0;

    // Parallel calculation
    start = clock();
    for (i = myid; i < n; i += numprocs) {
        y = c + j * dy + dy / 2.0;
        for (j = myid; j < m; j += numprocs) {
            x = a + i * dx + dx / 2.0;
            sum += f(x, y);
        }
        gsum = dx * sum;
    }
    gsum = dy * sum;

    // MPI Reduction
    MPI_Reduce(&sum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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

        // Compute the integral result
        result *= dx * dy;

        printf("Integral result = %.15f\n", result);
    }

    // MPI Finalization
    MPI_Finalize();

    return 0;
}