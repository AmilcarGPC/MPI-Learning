/*
PROGRAM: Calculate Heat Equation                                                        
HOW TO RUN :                                            
$ mpicxx Project1_MPI.cpp -o run_Project1_MPI
$ mpiexec -n 4 ./run_Project1_MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <cstdio>
#define pi 3.14159265359
#define PARAVIEW_MAX_GRID 150

void save_paraview(int Nx, int Ny, double *x, double *y, double t, double **Unew, double **PG, int nrec);

int main(int argc, char *argv[]) {
    int i, j, n, index;
    int NxG, NyG, Nt, Nx, Ny, NN, NNF, nrec, iIni, iFin;
    double kx, ky, xI, xF, yI, yF, tF, tI, Dx, Dy, Dt, t;
    double rx, ry, aC, aN, aS, aE, aW;
    double suma, sumaglob;
    double inicio, final, tiempo;
    double *x, *y, *xG, *yG;

    double **Uold = nullptr;
    double **Unew = nullptr;
    double **UG = nullptr;
    double **PG = nullptr;
    
    int ierr, numtasks, taskid, tipo_col, tipo_block, etiqueta;
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
    NxG = 150;   // Numero de puntos en x (GLOBAL)
    NyG = 150;   // Numero de puntos en y (GLOBAL)

    // Discretización
    Dx = (xF - xI) / (NxG - 1);
    Dy = (yF - yI) / (NyG - 1);
    Dt = (tF - tI) / (Nt - 1);
    rx = kx * (Dt / (Dx * Dx));
    ry = ky * (Dt / (Dy * Dy));

    xG = (double *)malloc(NxG * sizeof(double));
    yG = (double *)malloc(NyG * sizeof(double));
    UG = (double **)malloc(NxG * sizeof(double *));
    PG = (double **)malloc(NxG * sizeof(double *));
    for (i = 0; i < NxG; i++){
        UG[i] = (double *)malloc(NyG * sizeof(double));
        PG[i] = (double *)malloc(NyG * sizeof(double));
    }

    for (i = 0; i < NxG; i++)
        for (j = 0; j < NyG; j++)
            PG[i][j] = 0.0;

    // Malla
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

    // MPI: Topologia (cartesiana) 
    vecino[0] = MPI_PROC_NULL;
    vecino[1] = MPI_PROC_NULL;
    if (taskid > 0)
        vecino[0] = taskid - 1;
    if (taskid < numtasks - 1)
        vecino[1] = taskid + 1;

    // MPI: División del dominio

    // División del dominio en x
    NN = floor(1.0 * NxG / numtasks);
    NNF = NxG - NN * taskid;
    Nx = NN + 2;
    if (taskid == 0)
        Nx = NN + 1;
    if (taskid == numtasks - 1)
        Nx = NNF + 1;
    // División del dominio en y
    Ny = NyG;

    // MPI: Indices locales y globales
    index_global = (int *)malloc(Nx * sizeof(int));

    if (taskid == 0) {
        for (i = 0; i < Nx; i++)
            index_global[i] = i;
    } else {
        for (i = 0; i < Nx; i++)
            index_global[i] = taskid * NN + i - 1;
    }

    // Imprimir
    if (taskid == 0) {
        printf("\n");
        printf("(NxG,NyG)=%d %d\n", NxG, NyG);
        printf("Criterio CFL < 1/2 : %.15f\n", rx + ry);
        printf("\n");
    }

    printf("taskid = %d, NN = %d, Nx_local = %d, index global = %d %d\n", taskid, NN, Nx, index_global[0], index_global[Nx - 1]);

    // MPI: Variables locales
    x = (double *)malloc(Nx * sizeof(double));
    y = (double *)malloc(Ny * sizeof(double));
    Uold = (double **)malloc(Nx * sizeof(double *));
    Unew = (double **)malloc(Nx * sizeof(double *));
    for (i = 0; i < Nx; i++) {
        Uold[i] = (double *)malloc(Ny * sizeof(double));
        Unew[i] = (double *)malloc(Ny * sizeof(double));
    }

    // Malla
    for (i = 0; i < Nx; i++) {
        index = index_global[i];
        x[i] = xG[index];
    }
    for (j = 0; j < Ny; j++) {
        y[j] = yG[j];
    }

    // Condiciones Iniciales
    for (i = 0; i < Nx; i++) {
        index = index_global[i];
        for (j = 0; j < Ny; j++) {
            Uold[i][j] = UG[index][j];
        }
    }

    // Actualizacion
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            Unew[i][j] = Uold[i][j];
        }
    }

    // <----------< PROGRAMA PRINCIPAL >---------->
    //Guardar condiciones iniciales y de frontera
    nrec = 0;
    if (taskid == 0){
        save_paraview(NxG,NyG,xG,yG,t,UG,PG,nrec);
    }

    // MPI: Tipo de vectores (comunicaciones)
    MPI_Type_vector(1, Ny, 0, MPI_DOUBLE_PRECISION, &tipo_col);
    MPI_Type_vector(NN*Ny,1,1, MPI_DOUBLE_PRECISION, &tipo_block);
    MPI_Type_commit(&tipo_col);
    MPI_Type_commit(&tipo_block);

    // MPI: Loop de calculos
    inicio = MPI_Wtime();

    aC = 1 - 2.0 * (rx + ry);
    aE = rx;
    aW = rx;
    aS = ry;
    aN = ry;
    for (n = 0; n < Nt; n++) {
        // Nuevos valores
        for (i = 1; i < Nx - 1; i++) {
            for (j = 1; j < Ny - 1; j++) {
                Unew[i][j] = aC * Uold[i][j] +
                             aW * Uold[i - 1][j] +
                             aE * Uold[i + 1][j] +
                             aS * Uold[i][j - 1] +
                             aN * Uold[i][j + 1];
            }
        }

        // Comunicacion
        MPI_Sendrecv(&Unew[1][0], 1, tipo_col, vecino[0], etiqueta,
                         &Unew[0][0], 1, tipo_col, vecino[0], etiqueta,
                         MPI_COMM_WORLD, &statut);
        MPI_Sendrecv(&Unew[Nx - 2][0], 1, tipo_col, vecino[1], etiqueta,
                         &Unew[Nx - 1][0], 1, tipo_col, vecino[1], etiqueta,
                         MPI_COMM_WORLD, &statut);

        // Actualización
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                Uold[i][j] = Unew[i][j];
            }
        }

        // Guardar resultados cada 1000 iteraciones
        if (n % 1000 == 0) {
            t = tI + n * Dt;
            nrec++;
            
            if (Nx <= PARAVIEW_MAX_GRID) {
                if (taskid < numtasks - 1) {
                    iIni = 1;
                    if (taskid == 0)
                        iIni = 0;
                    iFin = iIni + NN - 1;
                    for (int c = 0; c < NN; c++)
                        MPI_Send(&Unew[iIni+c][0], 1, tipo_col, numtasks - 1, etiqueta, MPI_COMM_WORLD);
                }

                MPI_Barrier(MPI_COMM_WORLD);
                if (taskid == numtasks-1) {
                    double received_column[Ny];

                    iIni = taskid * NN - 1;
                    iFin = iIni + Nx;
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 1; i < Nx; i++) {
                            UG[iIni + i][j] = Unew[i][j];
                        }
                    }
                    for (int i = 0; i < numtasks-1; i++) {
                        iIni = i * NN;
                        iFin = iIni + NN - 2;
                        for (int c = 0; c < NN; c++){
                            MPI_Recv(&received_column[0], 1, tipo_block, i, etiqueta, MPI_COMM_WORLD, &statut);
                            for (int j = 0; j < Ny; j++) {
                                //printf("taskid: %d, column[%d][%d]: %f\n",i,c,j,received_column[j]);
                                UG[iIni + c][j] = received_column[j];
                            }
                        }
                            
                        for (int j = 0; j < Ny; j++) {
                            for (int k = 0; k < iFin - iIni + 1; k++) {
                                PG[iIni + k][j] = i;
                            }
                        }
                    }
                    save_paraview(NxG, NyG, xG, yG, t, UG, PG, nrec);
                }

            }
        }
    }
    printf("El procesador %d terminó su chamba\n",taskid);
    final = MPI_Wtime();

    // <----------< VERIFICAR RESULTADOS >---------->
    // Imprimir suma de los elementos para verificar
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

    MPI_Allreduce(&suma, &sumaglob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    suma = sumaglob;

    printf("taskid=%d, Suma total = %f\n", taskid, suma);

    // Tiempo final
    tiempo = final - inicio;
    printf("tiempo=%f, taskid=%d\n", tiempo, taskid);

    // <----------< FINALIZACIÓN >---------->
    MPI_Type_free(&tipo_col);
    MPI_Type_free(&tipo_block);
    MPI_Finalize();

    // Liberar memoria
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
    for (i = 0; i < NxG; i++)
        free(PG[i]);
    free(PG);
    
    return 0;
}

void save_paraview(int Nx, int Ny, double *x, double *y, double t, double **Unew, double **PG, int nrec) {
    int i, j;
    FILE *file;
    char filen[100]; // Aumenta el tamaño del nombre del archivo para incluir la ruta completa

    // Construcción de la ruta completa del archivo
    sprintf(filen, "%s/P-%04d.vtk", "paraview", nrec);

    // Apertura del archivo de salida
    file = fopen(filen, "w");

    // Escritura de metadatos en el archivo
    fprintf(file, "# vtk DataFile Version 2.0\n");
    fprintf(file, "Sample rectilinear grid\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET RECTILINEAR_GRID\n");
    fprintf(file, "DIMENSIONS %d %d 1\n", Nx, Ny);

    // Escritura de coordenadas en el archivo
    fprintf(file, "X_COORDINATES %d float\n", Nx);
    for (i = 0; i < Nx; i++) {
        fprintf(file, "%e\n", x[i]);
    }

    fprintf(file, "Y_COORDINATES %d float\n", Ny);
    for (j = 0; j < Ny; j++) {
        fprintf(file, "%e\n", y[j]);
    }

    fprintf(file, "Z_COORDINATES 1 float\n");
    fprintf(file, "0.0\n");

    // Escritura del número de puntos
    fprintf(file, "POINT_DATA %d\n", Nx * Ny);

    // Escritura de los datos del campo en el archivo
    fprintf(file, "SCALARS U float\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            fprintf(file, "%e\n", Unew[i][j]);
        }
    }

    fprintf(file, "SCALARS Proc float\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            fprintf(file, "%e ", PG[i][j]);
        }
        fprintf(file, "\n");
    }

    // Cierre del archivo de salida
    fclose(file);

    // Impresión de información en la consola
    printf("\nParaview No. %d, time: %.3f\n", nrec, t);
}
