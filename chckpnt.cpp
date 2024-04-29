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
    int divx = 3;
    int divy = 3;
    int i, j, n, index,i_index,j_index;
    int NxG, NyG, Nt, Nx, Ny, NN, NNF, nrec, iIni, iFin;
    double kx, ky, xI, xF, yI, yF, tF, tI, Dx, Dy, Dt, t;
    double rx, ry, aC, aN, aS, aE, aW;
    double suma, sumaglob;
    double inicio, final, tiempo;
    double *x, *y, *xG, *yG;

    int Ndom,NdomF;
    int Mdom,MdomF;
    int jIni,jFin;

    double **Uold = nullptr;
    double **Unew = nullptr;
    double **UG = nullptr;
    double **PG = nullptr;
    
    int ierr, numtasks, taskid, tipo_col, tipo_block, etiqueta;
    int CORDS[2] = {MPI_PROC_NULL, MPI_PROC_NULL};
    int vecino[4];
    int *index_global_i;
    int *index_global_j;

    int **t_index_global_i;
    int **t_index_global_j;

    const int ndims = 2;
    int dims_vec[ndims]; 
    int periodicite[ndims]; 
    int reorganisation;
    int coords[ndims];
    int pcoords[ndims];
    MPI_Comm comm2D;   

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

    Nt = 100000; // Numero de pasos en t
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
        for (j = 0; j < NyG; j++){
            UG[i][j] = 3.0 * sin(pi * xG[i] + pi * yG[j]) * sin(pi * xG[i] + pi * yG[j]);
        }
            
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

    if (divx*divy==numtasks){
        if (taskid == 0)
            printf("\n CORRECTO: numtasks = divx*divy = %d \n",divx*divy);
        }
    else {
        if (taskid == 0)
            printf("\n INCORRECTO: numtasks = %d  y  divx*divy = %d \n",
        numtasks,divx*divy);
        return 0;
    } 

    /*if (taskid == 0){
        printf("Matriz:\n");
        for (int i = 0; i < NxG; i++) {
            for (int j = 0; j < NyG; j++) {
                printf("%4f ", UG[i][j]); // Ajuste el ancho de campo según sea necesario
            }
            printf("\n");
        }
    }*/

    // MPI: Topologia (cartesiana)
    dims_vec[0] = divx;                                                   
    dims_vec[1] = divy;                                               
    reorganisation = 0;                                         
    periodicite[0] = 0; 
    periodicite[1] = 0;

    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims_vec,              
                    periodicite,reorganisation,
                    &comm2D); 
                        
    //==========================================
    // MPI [3]: Coordenadas de la malla topologica
    MPI_Cart_coords(comm2D,taskid,ndims,coords);
        // MPI: División del dominio
    //-----------------
    Ndom  = (int)(1.0*NxG/divx);
    NdomF = NxG - Ndom*(divx-1);
    //-----------------
    Nx = Ndom;
    if (coords[0]==divx-1){Nx = NdomF;}
    //-----------------
    iIni = coords[0]*Ndom;
    iFin = iIni + Nx;
    //==========================================

    //==========================================
    // MPI [3]: Division del trabajo  
    //-----------------
    Mdom  = (int)(1.0*NyG/divy);
    MdomF = NyG - Mdom*(divy-1);
    //-----------------
    Ny = Mdom ;
    if (coords[1] == divy-1){
        Ny = MdomF;
    }

    //-----------------
    jIni = coords[1]*Mdom;
    jFin = jIni + Ny; 

    int destino;

    vecino[0] = MPI_PROC_NULL;
    vecino[1] = MPI_PROC_NULL;
    vecino[2] = MPI_PROC_NULL;
    vecino[3] = MPI_PROC_NULL;
    
    if (coords[0] - 1 >= 0) {
        CORDS[0] = coords[0] - 1;
        CORDS[1] = coords[1];
        Nx++;
        MPI_Cart_rank(comm2D,CORDS, &destino);
        vecino[0] = destino;
        //printf("dest1: %d\n",destino);
    }
    
    if (coords[0] + 1 < divx) {
        CORDS[0] = coords[0] + 1;
        CORDS[1] = coords[1];
        Nx++;
        MPI_Cart_rank(comm2D,CORDS,&destino);
        vecino[1] = destino;
        //printf("dest2: %d\n",destino);
    }
    if (coords[1] - 1 >= 0) {
        CORDS[0] = coords[0];
        CORDS[1] = coords[1] - 1;
        Ny++;
        MPI_Cart_rank(comm2D,CORDS, &destino);
        vecino[2] = destino;
        //printf("dest3: %d\n",destino);
    }
    if (coords[1] + 1 < divy) {
        Ny++;
        CORDS[0] = coords[0];
        CORDS[1] = coords[1] + 1;
        MPI_Cart_rank(comm2D,CORDS, &destino);
        vecino[3] = destino;
        //printf("dest4: %d\n",destino);
    }
    //printf("B| taskid %d, Ny: %d\n",taskid, Ny);
    // MPI: Indices locales y globales
    index_global_i = (int *)malloc(Nx * sizeof(int));
    index_global_j = (int *)malloc(Ny * sizeof(int));

    t_index_global_i = (int **)malloc(numtasks * sizeof(int *));
    for (i = 0; i < numtasks; i++) {
        t_index_global_i[i] = (int *)malloc(Nx * sizeof(int));
    }

    t_index_global_j = (int **)malloc(numtasks* sizeof(int *));
    for (i = 0; i < numtasks; i++) {
        t_index_global_j[i] = (int *)malloc(Ny * sizeof(int));
    }

    for (i = 0; i < numtasks; i++) {
        MPI_Cart_coords(comm2D,i,ndims,pcoords);
        //printf("A| tasdkid: %d. coordenadas (%d,%d).\n",i,pcoords[0],pcoords[1]);
        if (pcoords[0] == 0) {
            for (j = 0; j < Nx; j++)
                t_index_global_i[i][j] = j;
        } else {
            for (j = 1; j < Nx + 1; j++)
                t_index_global_i[i][j-1] = pcoords[0] * Ndom + j - 1;
        }
    }

    for (i = 0; i < numtasks; i++) {
        MPI_Cart_coords(comm2D,i,ndims,pcoords);
        //printf("A| tasdkid: %d. coordenadas (%d,%d).\n",i,pcoords[0],pcoords[1]);
        if (pcoords[1] == 0) {
            for (j = 0; j < Ny; j++)
                t_index_global_j[i][j] = j;
        } else {
            for (j = 1; j < Ny + 1; j++)
                t_index_global_j[i][j-1] = pcoords[1] * Mdom + j - 1;
        }
    }

    //printf("tsk: %d, [%d][%d]\n",taskid, t_index_global_i[taskid][0], t_index_global_j[taskid][0]);

    if (coords[0] == 0) {
        for (i = 0; i < Nx; i++)
            index_global_i[i] = i;
    } else {
        for (i = 0; i < Nx; i++)
            index_global_i[i] = coords[0] * Ndom + i - 1;
    }

    if (coords[1] == 0) {
        for (j = 0; j < Ny; j++)
            index_global_j[j] = j;
    } else {
        for (j = 0; j < Ny; j++)
            index_global_j[j] = coords[1] * Mdom + j - 1;
    }

    // Imprimir
    if (taskid == 0) {
        printf("\n");
        printf("(NxG,NyG)=%d %d\n", NxG, NyG);
        printf("Criterio CFL < 1/2 : %.15f\n", rx + ry);
        printf("\n");
    }
    printf("taskid = %d, Ndom = %d, Mdom = %d, Nx_local = %d, Ny_local = %d, index global i = [%d->%d], index global j = [%d->%d]\n", taskid, Ndom, Mdom, Nx, Ny, index_global_i[0], index_global_i[Nx - 1], index_global_j[0], index_global_j[Ny - 1]);
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
        i_index = index_global_i[i];
        x[i] = xG[i_index];
    }
    for (j = 0; j < Ny; j++) {
        j_index = index_global_j[j];
        y[j] = yG[j_index];
    }
    //printf("A| tasdkid: %d\n",taskid);
    // Condiciones Iniciales
        for (i = 0; i < Nx; i++) {
            i_index = index_global_i[i];
            for (j = 0; j < Ny; j++) {
                j_index = index_global_j[j];
                //printf("UG[%d][%d] = %f\n",i_index,j_index,UG[i_index][j_index]);
                Uold[i][j] = UG[i_index][j_index];
            }
        }
    //printf("B| tasdkid: %d\n",taskid);
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
    MPI_Type_vector(1, MdomF, 0, MPI_DOUBLE_PRECISION, &tipo_col);
    MPI_Type_vector(NN*Ny,1,1, MPI_DOUBLE_PRECISION, &tipo_block);
    MPI_Type_commit(&tipo_col);
    MPI_Type_commit(&tipo_block);
    MPI_Type_commit(&tipo_col);

    // MPI: Loop de calculos
    inicio = MPI_Wtime();

    aC = 1 - 2.0 * (rx + ry);
    aE = rx;
    aW = rx;
    aS = ry;
    aN = ry;
    //printf("C| tasdkid: %d\n",taskid);

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

        // North exchange
        if (vecino[3] != MPI_PROC_NULL) {
            MPI_Request req[2];
            for (int i = 0; i < Nx; i++) {
                MPI_Isend(&Unew[i][Ny - 2], 1, MPI_DOUBLE_PRECISION, vecino[3], etiqueta, comm2D, &req[0]);
                MPI_Irecv(&Unew[i][Ny - 1], 1, MPI_DOUBLE_PRECISION, vecino[3], etiqueta, comm2D, &req[1]);
                MPI_Waitall(2, req, &statut);
            }
        }
        if (vecino[2] != MPI_PROC_NULL) {
            MPI_Request req1[2];
            for (int i = 0; i < Nx; i++) {
                MPI_Isend(&Unew[i][1], 1, MPI_DOUBLE_PRECISION, vecino[2], etiqueta, comm2D, &req1[0]);
                MPI_Irecv(&Unew[i][0], 1, MPI_DOUBLE_PRECISION, vecino[2], etiqueta, comm2D, &req1[1]);
                MPI_Waitall(2, req1, &statut);
            }
        }

        MPI_Sendrecv(&Unew[1][0], 1, tipo_col, vecino[0], etiqueta,
                         &Unew[0][0], 1, tipo_col, vecino[0], etiqueta,
                         comm2D, &statut);
        MPI_Sendrecv(&Unew[Nx - 2][0], 1, tipo_col, vecino[1], etiqueta,
                         &Unew[Nx - 1][0], 1, tipo_col, vecino[1], etiqueta,
                         comm2D, &statut);

        //printf("C| tasdkid: %d\n",taskid);
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
            //if (taskid == 0)
                //printf("%d tiempo: %f\n",nrec,t);

            if (Nx <= PARAVIEW_MAX_GRID) {
                if (taskid < numtasks-1) {
                    int lNx;
                    lNx = Ndom;
                    if (coords[0]==divx-1){
                        lNx = NdomF;
                    }
                    iIni = 1;
                    if (coords[0] == 0)
                        iIni = 0;

                    for (int c = 0; c < lNx; c++){
                        MPI_Send(&Unew[iIni+c][0], 1, tipo_col, numtasks - 1, etiqueta, MPI_COMM_WORLD);
                    }             
                }

                MPI_Barrier(comm2D);
                if (taskid == numtasks-1) {
                    double received_column[Ny];
                    iIni = t_index_global_i[taskid][0];
                    jIni = t_index_global_j[taskid][0];
                    
                    if (coords[0] == 0){
                        for (int j = 0; j < MdomF; j++) {
                            for (int i = 0; i < NdomF; i++) {
                                UG[iIni+i][jIni+j-1] = Unew[i][j];
                            }
                        }
                    }
                    
                    if (coords[1] == 0){
                        for (int j = 0; j < MdomF; j++) {
                            for (int i = 0; i < NdomF; i++) {
                                UG[iIni+i][jIni+j] = Unew[i+1][j];
                            }
                        }
                    }

                    if ((coords[0] != 0) && (coords[1] != 0)){
                        for (int j = 1; j < MdomF + 1; j++) {
                            for (int i = 1; i < NdomF + 1; i++) {
                                UG[iIni+i-1][jIni+j-2] = Unew[i][j-1];
                            }
                        }
                    }

                    for (int i = 0; i < numtasks-1; i++) {
                        int lNx;
                        int lNy;
                        int lcoords[ndims];

                        MPI_Cart_coords(comm2D,i,ndims,lcoords);
                        lNx = Ndom;
                        if (lcoords[0] == divx-1){lNx = NdomF;}
                        lNy = Mdom ;
                        if (lcoords[1] == divy-1){lNy = MdomF;}
                        iIni = t_index_global_i[i][0];
                        jIni = t_index_global_j[i][0];
                        iFin = t_index_global_i[i][lNx-1];
                        jFin = t_index_global_j[i][lNy-1];

                        for (int c = 0; c < lNx; c++){
                            MPI_Recv(&received_column[0], 1, tipo_col, i, etiqueta, MPI_COMM_WORLD, &statut);
                            
                            if (lcoords[1] == 0){
                                for (int j = 0; j < lNy; j++) {
                                    UG[iIni + c][jIni + j] = received_column[j];
                                }
                            } else {
                                for (int j = 1; j < lNy + 1; j++) { //lNy + 1 causa desbordes
                                    UG[iIni + c][jIni + j - 2] = received_column[j-1];
                                }
                            }
                        }

                        for (int j = 0; j < Ny; j++) {
                            for (int k = 0; k < iFin - iIni + 1; k++) {
                                PG[iIni + k][j] = i;
                            }
                        }
                    }

                    /*if (nrec == 2)
                        for (int h = 0; h < NxG; h++)
                            printf("UG[%d][%d] = %f\n",h,100,UG[h][149]);
*/
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
    free(index_global_i);
    free(index_global_j);
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
            //fprintf(file, "Unew[%d][%d]: %e\n", i,j,Unew[i][j]);
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
