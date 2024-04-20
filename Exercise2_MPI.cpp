/*
PROGRAM: Sum of elements of a vector                                                            !
HOW TO RUN :                                            
$ mpicxx Exercise2_MPI.cpp -o Exercise2_MPI          
$ mpiexec -n 4 ./Exercise2_MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

double *generateRandomArray(int numberOfElements);
double *generateSinArray(int numberOfElements);

int main()
{
   int N = 100000005;
   int i,iIni,iFin;
   double suma,cpu_time_used;
   double *x;      
   //==================
   // MPI [1]: definiciones
   int numtasks, taskid;
   int ndims = 2;
   int dims[ndims];
   int periodicite[ndims];
   int reorganisation;
   int comm1D;
   int Nx,Ndom;
   double sumaglob;
   double sumaglobito;
   double *sumaux;
   //==================

   //__________________________________________
   // Vector de N elementos
   //x = generateRandomArray(N); /* Opcion 1: Arbitrario */
   x = generateSinArray(N);      /* Opcion 2: Integral del seno */
      
   //==========================================
   // MPI [2]: inicio
   //-----------------------
   MPI_Init(NULL,NULL);                     /* Initializacion*/
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks); /* Numero de procesos */
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);   /* Identificador de cada proceso */
   //==========================================

   iIni = 0;
   iFin = N;
   //==========================================
   // MPI [3]: Division del trabajo
   Ndom = (int)(1.0*N/numtasks);
   if (numtasks == 1) {
       Nx = N;
   }
   else {
       if(taskid==numtasks-1){
           Nx = N - Ndom*(numtasks-1);
       }
       else {
           Nx = Ndom;
       }
   }
   iIni = taskid*Ndom;
   iFin = iIni + Nx;
   //==========================================
    
   //__________________________________________
   // Operaciones: sumar sus componentes
   clock_t start, end;
   start = clock();
   //----------------------- 
   suma = 0.0;
   for (i = iIni; i < iFin; i++) {
        suma = suma + x[i];
   }

   //==========================================
   // MPI    
   sumaux = (double *)malloc(numtasks*sizeof(double));
   //==========================================

   //==========================================
   // MPI [4]: sumar componentes de cada procesador
   if (numtasks > 1){
      if (taskid != 0) {
            MPI_Send(&suma,1,MPI_DOUBLE_PRECISION,0,1000,MPI_COMM_WORLD);
        } else {
            printf(" \n");
            sumaglob = suma; 
            for (int i = 1; i < numtasks; i++){
                MPI_Recv(&sumaglobito,1,MPI_DOUBLE_PRECISION,i,1000,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("Process 0 received number %f from process %d\n",sumaglobito,i);
                sumaglob += sumaglobito;
            }
            printf("\nTotal sum: %f \n\n",sumaglob);
        }
   } else {
       sumaglob = suma;
   }
   //==========================================

   //-----------------------
   end = clock();
   cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC;
   
   //==========================================
   // MPI [5]: sumar componentes de cada procesador
   // MPI_Allreduce(&suma,&sumaglob,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
   //==========================================

   //__________________________________________
   // Imprimir
   printf("Procesor %i de %i: \n \
           Nx    = %i  \n \
           iIni  = %i  \n \
           iFin  = %i  \n \
           suma  = %f  \n \
           sumaT = %f  \n \
           time  = %lf \n \n",taskid,numtasks,Nx,iIni,iFin,suma,sumaglob,cpu_time_used);

   //__________________________________________
   // Liberar memoria
   free(x);
   free(sumaux);

   //==========================================
   // MPI [5]: final
   MPI_Finalize();
   //==========================================

   return 0;
}

//sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
// Vector de numeros aleatorios entre 0 y 100 

double *generateRandomArray(int N){
   double *x;
   srand(time(NULL));
   x = (double *)malloc(N*sizeof(double));
   for (int i = 0; i < N; i++){
       x[i] = rand() % 100 + 1;
       //printf("%f\n",x[i]);
   }
   return x;
}

//sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
// Integral del seno
 
double *generateSinArray(int N){
   double *x;
   double xx,Dx;
   x = (double *)malloc(N*sizeof(double));
   Dx = 1.0/(N-1);   
   for (int i = 0; i < N; i++){
        xx = i*Dx;
        x[i] = sin(xx)/N;
        //printf("%f\n",x[i]);
   }
   return x;
}

//sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss