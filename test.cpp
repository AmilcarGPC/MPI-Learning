/*
PROGRAM: Bidimentional Calculate Heat Equation                                                        
HOW TO RUN :                                            
$ mpicxx test.cpp -o run_test
$ mpiexec -n 9 ./run_test 3 3
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <cstdio>
#define pi 3.14159265359
#define PARAVIEW_MAX_GRID 150

void save_paraview(int Nx, int Ny, double *x, double *y, double t, double *Unew, double *PG, int nrec);

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Uso: %s <numero1> <numero2>\n", argv[0]);
        return 1;
    }

    // Convertir los argumentos de cadena a números enteros
    int divx = atoi(argv[1]);
    int divy = atoi(argv[2]);
    int i, j, n, i_index, j_index, up;
    int NxG, NyG, Nt, Nx, Ny, nrec;
    double kx, ky, xI, xF, yI, yF, tF, tI, Dx, Dy, Dt, t;
    double rx, ry, aC, aN, aS, aE, aW;
    double suma, sumaglob;
    double inicio, final, tiempo, tiempoglob;
    double *x, *y, *xG, *yG;

    double avg5times[5];
    double seqTime = 520.019127;

    int Ndom, NdomF;
    int Mdom, MdomF;
    int iIni, iFin, jIni, jFin;
    int GiIni, GjIni;

    double *Uold = nullptr;
    double *Unew = nullptr;
    double *UG = nullptr;
    double *PG = nullptr;
    
    int lNx, lNy;
    int numtasks, taskid, tipo_block, tipo_row, tipo_col, etiqueta;
    int tipo_block_rcv_1, tipo_block_rcv_2, tipo_block_rcv_3;

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
    UG = (double *)malloc(NxG * NyG * sizeof(double));
    PG = (double *)malloc(NxG * NyG * sizeof(double));
    up = 0;
    for (i = 0; i < NxG; i++)
        for (j = 0; j < NyG; j++)
            PG[up++] = 0.0;

    // Malla global
    for (i = 0; i < NxG; i++)
        xG[i] = xI + i * Dx;
    for (j = 0; j < NyG; j++)
        yG[j] = yI + j * Dy;

    // Condiciones Iniciales
    up = 0;
    for (i = 0; i < NxG; i++)
        for (j = 0; j < NyG; j++){
            UG[up++] = 3.0 * sin(pi * xG[i] + pi * yG[j]) * sin(pi * xG[i] + pi * yG[j]);
        }
            
    // Condiciones de frontera
    for (j = 0; j < NyG; j++) {
        UG[j] = 2.0;     // Oeste
        UG[(NxG - 1)*NyG+j] = 1.0; // Este
    }
    for (i = 0; i < NyG; i++) {
        UG[i*NyG] = 1.0;     // Sur
        UG[i*NyG+NyG - 1] = 3.0; // Norte
    }
    // <----------< COMPONENTES MPI >---------->
    // MPI: Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    if (taskid == 0)
        printf("\nConfig: %d × %d: %d",divx, divy, divx*divy);
    for (int ephocs = 0; ephocs < 5; ephocs++){
        MPI_Status statut;

        if (divx*divy!=numtasks){
            MPI_Barrier(MPI_COMM_WORLD);
            return 0;
        } 

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

        // MPI: División del dominio en X
        Ndom  = NxG / divx;
        NdomF = NxG - Ndom * (divx - 1);
        Nx = (coords[0] < divx - 1) ? Ndom : NdomF;

        // MPI: División del dominio en Y
        Mdom  = NyG / divy;
        MdomF = NyG - Mdom * (divy - 1);
        Ny = (coords[1] < divy - 1) ? Mdom : MdomF;


        // MPI: Encontrar vecinos del proceso
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
        }
        
        if (coords[0] + 1 < divx) {
            CORDS[0] = coords[0] + 1;
            CORDS[1] = coords[1];
            Nx++;
            MPI_Cart_rank(comm2D,CORDS,&destino);
            vecino[1] = destino;
        }

        if (coords[1] - 1 >= 0) {
            CORDS[0] = coords[0];
            CORDS[1] = coords[1] - 1;
            Ny++;
            MPI_Cart_rank(comm2D,CORDS, &destino);
            vecino[2] = destino;
        }

        if (coords[1] + 1 < divy) {
            Ny++;
            CORDS[0] = coords[0];
            CORDS[1] = coords[1] + 1;
            MPI_Cart_rank(comm2D,CORDS, &destino);
            vecino[3] = destino;
        }

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
            if (pcoords[1] == 0) {
                for (j = 0; j < Ny; j++)
                    t_index_global_j[i][j] = j;
            } else {
                for (j = 1; j < Ny + 1; j++)
                    t_index_global_j[i][j-1] = pcoords[1] * Mdom + j - 1;
            }
        }

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

        // MPI: Variables locales
        x = (double *)malloc(Nx * sizeof(double));
        y = (double *)malloc(Ny * sizeof(double));
        Uold = (double *)malloc(Nx * Ny * sizeof(double));
        Unew = (double *)malloc(Nx * Ny * sizeof(double));
        
        // MPI: Malla local
        for (i = 0; i < Nx; i++) {
            i_index = index_global_i[i];
            x[i] = xG[i_index];
        }
        for (j = 0; j < Ny; j++) {
            j_index = index_global_j[j];
            y[j] = yG[j_index];
        }
        
        // MPI: Condicion inicial local
        up = 0;
        for (i = 0; i < Nx; i++) {
            i_index = index_global_i[i];
            for (j = 0; j < Ny; j++) {
                j_index = index_global_j[j];
                Uold[up++] = UG[i_index*NxG+j_index];
            }
        }
            
        // MPI: Actualizacion
        up = 0;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                Unew[up] = Uold[up];
                up++;
            }
        }
        
        // <----------< PROGRAMA PRINCIPAL >---------->
        //MPI: Guardar condiciones iniciales y de frontera
        nrec = 0;
        
        // Inicio y fin en X y Y sin contar la extensión de dominio
        iIni = 1;
        if (coords[0] == 0)
            iIni = 0;
        jIni = 1;
        if (coords[1] == 0)
            jIni = 0;
        iFin = iIni+Ndom;
        if (coords[0] == divx-1)
            iFin = iIni+NdomF;
        jFin = jIni+Mdom;
        if (coords[1] == divy-1)
            jFin = jIni+MdomF;

        // MPI: Tipo de vectores (comunicaciones)
        MPI_Type_vector(1, Ny, 0, MPI_DOUBLE_PRECISION, &tipo_row);
        MPI_Type_vector(Nx, 1, Ny, MPI_DOUBLE_PRECISION, &tipo_col);
        MPI_Type_vector((iFin-iIni),(jFin-jIni),Ny,MPI_DOUBLE_PRECISION, &tipo_block);
        MPI_Type_commit(&tipo_block);
        MPI_Type_commit(&tipo_row);
        MPI_Type_commit(&tipo_col);

        // Tipo de vectores para recibir los bloques (hay 3 tipos)
        if (taskid == numtasks - 1){
            MPI_Type_vector(Ndom,Mdom,Mdom,MPI_DOUBLE_PRECISION, &tipo_block_rcv_1);
            MPI_Type_vector(NdomF,Mdom,Mdom,MPI_DOUBLE_PRECISION, &tipo_block_rcv_2);
            MPI_Type_vector(Ndom,MdomF,MdomF,MPI_DOUBLE_PRECISION, &tipo_block_rcv_3);
            MPI_Type_commit(&tipo_block_rcv_1);
            MPI_Type_commit(&tipo_block_rcv_2);
            MPI_Type_commit(&tipo_block_rcv_3);
        }

        // MPI: Loop de calculos
        inicio = MPI_Wtime();

        aC = 1 - 2.0 * (rx + ry);
        aE = rx;
        aW = rx;
        aS = ry;
        aN = ry;

        // Función iterativa para calcular el calor
        for (n = 0; n < Nt; n++) {
            // Nuevos valores
            for (i = 1; i < Nx - 1; i++) {
                for (j = 1; j < Ny - 1; j++) {
                    Unew[i*Ny+j] = aC * Uold[i*Ny+j] +
                                aW * Uold[(i-1)*Ny+j] +
                                aE * Uold[(i+1)*Ny+j] +
                                aS * Uold[i*Ny+j-1] +
                                aN * Uold[i*Ny+j+1];         
                }
            }

            // Comunicación
            MPI_Sendrecv(&Unew[Ny], 1, tipo_row, vecino[0], etiqueta,
                            &Unew[0], 1, tipo_row, vecino[0], etiqueta,
                            comm2D, &statut);
                            
            MPI_Sendrecv(&Unew[(Nx - 2)*Ny], 1, tipo_row, vecino[1], etiqueta,
                            &Unew[(Nx - 1)*Ny], 1, tipo_row, vecino[1], etiqueta,
                            comm2D, &statut);
                            
            MPI_Sendrecv(&Unew[1], 1, tipo_col, vecino[2], etiqueta,
                            &Unew[0], 1, tipo_col, vecino[2], etiqueta,
                            comm2D, &statut);
            MPI_Sendrecv(&Unew[Ny - 2], 1, tipo_col, vecino[3], etiqueta,
                            &Unew[Ny - 1], 1, tipo_col, vecino[3], etiqueta,
                            comm2D, &statut);
            
            // Actualización
            up = 0;
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    Uold[up] = Unew[up];
                    up++;
                }
            }
        }
        final = MPI_Wtime();
        // <----------< VERIFICAR RESULTADOS >---------->
        // Imprimir suma de los elementos para verificar
        suma = 0.0;
        for (i = iIni; i < iFin; i++) {
            for (j = jIni; j < jFin; j++) {
                suma += Unew[i*Ny+j];
            }
        }
        tiempo = final - inicio;
        MPI_Allreduce(&suma, &sumaglob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&tiempo, &tiempoglob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        suma = sumaglob;
        if (taskid == 0){
            printf("\nEpoch: %d | avgT: %f | Suma total = %f |",ephocs, tiempoglob/numtasks, suma);
            avg5times[ephocs] = tiempoglob/numtasks;
        }

        // <----------< FINALIZACIÓN >---------->
        if (taskid == numtasks - 1){
            MPI_Type_free(&tipo_block_rcv_1);
            MPI_Type_free(&tipo_block_rcv_2);
            MPI_Type_free(&tipo_block_rcv_3);
        }

        MPI_Type_free(&tipo_block);
        MPI_Type_free(&tipo_row);
        MPI_Type_free(&tipo_col);
    }
    MPI_Finalize();
    
    if (taskid == 0){
        tiempo = 0.0;
        for(int i = 0; i < 5; i++){
            tiempo += avg5times[i];
        }
        tiempo = tiempo/5;
        printf("\n>>> Statistics: avg: %f | speedUp: %f | efficiency: %f |",tiempo,seqTime/tiempo,(seqTime/tiempo)/(divx*divy));
    }
        

    // Liberar memoria
    free(x);
    free(y);
    free(xG);
    free(yG);
    free(index_global_i);
    free(index_global_j);

    for (i = 0; i < numtasks; i++) {
        free(t_index_global_i[i]);
        free(t_index_global_j[i]);
    }
    free(t_index_global_i);
    free(t_index_global_j);
    free(Uold);
    free(Unew);
    free(UG);
    free(PG);
    
    return 0;
}

void save_paraview(int Nx, int Ny, double *x, double *y, double t, double *Unew, double *PG, int nrec) {
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
            fprintf(file, "%e\n", Unew[i*Ny+j]);
        }
    }

    fprintf(file, "SCALARS Proc float\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            fprintf(file, "%e ", PG[i*Nx+j]);
        }
        fprintf(file, "\n");
    }

    // Cierre del archivo de salida
    fclose(file);

    // Impresión de información en la consola
    printf("\nParaview No. %d, time: %.3f\n", nrec, t);
}
