/*
PROGRAM: Unidimentional Calculate Heat Equation                                                     
HOW TO RUN :                                            
$ mpicxx Project1_UniD_MPI.cpp -o run_Project1_UniD_MPI
$ mpiexec -n 9 ./run_Project1_UniD_MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define pi 3.14159265359

int main(int argc, char *argv[]) {
    int i, j, n, index;
    int NxG, NyG, Nt, Nx, Ny, NN, NNF, nrec, iIni, iFin;
    double kx, ky, xI, xF, yI, yF, tF, tI, Dx, Dy, Dt, t;
    double rx, ry, aC, aN, aS, aE, aW;
    double suma, sumaglob;
    double inicio, final, tiempo;
    double *x, *y, *xG, *yG;
    double **Uold, **Unew, **UG;
    int ierr, numtasks, taskid, tipo_row, etiqueta;
    int vecino[2];
    int *index_global;

    etiqueta = 2001;

    // <----------< INICIALIZACIÓN >---------->
    // Parametros
    kx = 1.0; // Termino difusivo en x
    ky = 1.0; // Termino difusivo en y

    xI = 0.0; // Inicio del dominio en x
    xF = 1.0; // Fin del dominio en x
    yI = 0.0; // Inicio del dominio en y
    yF = 1.0; // Fin del dominio en y

    tI = 0.0; // Tiempo inicial
    tF = 0.2; // Tiempo final

    Nt = 200000; // Numero de pasos en t
    NxG = 500;   // Numero de puntos en x (GLOBAL)
    NyG = 500;   // Numero de puntos en y (GLOBAL)

    // Discretización
    Dx = (xF - xI) / (NxG - 1);
    Dy = (yF - yI) / (NyG - 1);
    Dt = (tF - tI) / (Nt - 1);
    rx = kx * (Dt / (Dx * Dx));
    ry = ky * (Dt / (Dy * Dy));

    xG = (double *)malloc(NxG * sizeof(double));
    yG = (double *)malloc(NyG * sizeof(double));
    UG = (double **)malloc(NxG * sizeof(double *));
    for (i = 0; i < NxG; i++)
        UG[i] = (double *)malloc(NyG * sizeof(double));

    // Malla global
    for (i = 0; i < NxG; i++)
        xG[i] = xI + i * Dx;
    for (j = 0; j < NyG; j++)
        yG[j] = yI + j * Dy;

    // Condiciones Iniciales
    for (i = 0; i < NxG; i++)
        for (j = 0; j < NyG; j++)
            UG[i][j] = 3.0 * sin(pi * xG[i] + pi * yG[j]) * sin(pi * xG[i] + pi * yG[j]);

    // Condiciones de frontera
    for (j = 0; j < NyG; j++) {
        UG[0][j] = 2.0;     // Oeste
        UG[NxG - 1][j] = 1.0; // Este
    }
    for (i = 0; i < NxG; i++) {
        UG[i][0] = 1.0;     // Sur
        UG[i][NyG - 1] = 3.0; // Norte
    }

    // <----------< COMPONENTES MPI >---------->
    // MPI: Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Status statut;

    // MPI: Topologia
    vecino[0] = MPI_PROC_NULL;
    vecino[1] = MPI_PROC_NULL;
    if (taskid > 0)
        vecino[0] = taskid - 1;
    if (taskid < numtasks - 1)
        vecino[1] = taskid + 1;

    // MPI: Subdominio
    NN = floor(1.0 * NxG / numtasks);
    NNF = NxG - NN * taskid;
    Nx = NN + 2;
    if (taskid == 0)
        Nx = NN + 1;
    if (taskid == numtasks - 1)
        Nx = NNF + 1;
    Ny = NyG;

    index_global = (int *)malloc(Nx * sizeof(int));

    // MPI: Relación de índices 
    if (taskid == 0) {
        for (i = 0; i < Nx; i++)
            index_global[i] = i;
    } else {
        for (i = 0; i < Nx; i++)
            index_global[i] = taskid * NN + i - 1;
    }

    if (taskid == 0) {
        printf("\n");
        printf("(NxG,NyG)=%d %d\n", NxG, NyG);
        printf("Criterio CFL < 1/2 : %.15f\n", rx + ry);
        printf("\n");
    }

    printf("taskid = %d, NN = %d, Nx_local = %d, index global = %d %d\n", taskid, NN, Nx, index_global[0], index_global[Nx - 1]);

    x = (double *)malloc(Nx * sizeof(double));
    y = (double *)malloc(Ny * sizeof(double));
    Uold = (double **)malloc(Nx * sizeof(double *));
    Unew = (double **)malloc(Nx * sizeof(double *));
    for (i = 0; i < Nx; i++) {
        Uold[i] = (double *)malloc(Ny * sizeof(double));
        Unew[i] = (double *)malloc(Ny * sizeof(double));
    }

    // MPI: Malla local
    for (i = 0; i < Nx; i++) {
        index = index_global[i];
        x[i] = xG[index];
    }
    for (j = 0; j < Ny; j++) {
        y[j] = yG[j];
    }

    // MPI: Condicion inicial local
    for (i = 0; i < Nx; i++) {
        index = index_global[i];
        for (j = 0; j < Ny; j++) {
            Uold[i][j] = UG[index][j];
        }
    }

    // MPI: Actualización
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            Unew[i][j] = Uold[i][j];
        }
    }

    // MPI: Definición de tipos vector
    MPI_Type_vector(1, Ny, 0, MPI_DOUBLE_PRECISION, &tipo_row);
    MPI_Type_commit(&tipo_row);

    // <----------< PROGRAMA PRINCIPAL >---------->
    // MPI: Loop de calculos
    aC = 1 - 2.0 * (rx + ry);
    aE = rx;
    aW = rx;
    aS = ry;
    aN = ry;
    
    // Función iterativa para calcular el calor
    inicio = MPI_Wtime();
    for (n = 0; n < Nt; n++) {
        // MPI: Nuevos valores
        for (i = 1; i < Nx - 1; i++) {
            for (j = 1; j < Ny - 1; j++) {
                Unew[i][j] = aC * Uold[i][j] +
                             aW * Uold[i - 1][j] +
                             aE * Uold[i + 1][j] +
                             aS * Uold[i][j - 1] +
                             aN * Uold[i][j + 1];
            }
        }

        // MPI: Comunicación con procesos vecinos
        MPI_Sendrecv(&Unew[1][0], 1, tipo_row, vecino[0], etiqueta,
                         &Unew[0][0], 1, tipo_row, vecino[0], etiqueta,
                         MPI_COMM_WORLD, &statut);
        MPI_Sendrecv(&Unew[Nx - 2][0], 1, tipo_row, vecino[1], etiqueta,
                         &Unew[Nx - 1][0], 1, tipo_row, vecino[1], etiqueta,
                         MPI_COMM_WORLD, &statut);

        // MPI: Actualización
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                Uold[i][j] = Unew[i][j];
            }
        }
    }
    final = MPI_Wtime();

    // <----------< VERIFICAR RESULTADOS >---------->
    // MPI: Calculo de suma local
    iIni = 1;
    iFin = Nx - 1;
    if (taskid == 0)
        iIni = 0;
    if (taskid == numtasks - 1)
        iFin = Nx;

    suma = 0.0;
    for (i = iIni; i < iFin; i++) {
        for (j = 0; j < Ny; j++) {
            suma += Unew[i][j];
        }
    }

    // MPI: Sumar todas las sumas locales
    MPI_Allreduce(&suma, &sumaglob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // MPI: Suma global
    suma = sumaglob;

    printf("taskid=%d, Suma total = %f\n", taskid, suma);

    tiempo = final - inicio;
    printf("tiempo=%f, taskid=%d\n", tiempo, taskid);

    MPI_Type_free(&tipo_row);

    // <----------< FINALIZACIÓN >---------->
    // Liberar memoria
    MPI_Finalize();
    free(x);
    free(y);
    free(index_global);
    for (i = 0; i < Nx; i++) {
        free(Uold[i]);
        free(Unew[i]);
    }
    free(Uold);
    free(Unew);
    free(xG);
    free(yG);
    for (i = 0; i < NxG; i++)
        free(UG[i]);
    free(UG);
    
    return 0;
}