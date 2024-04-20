/*
PROGRAM: Vector addition
HOW TO RUN :
$ mpicxx Exercise1_MPI.cpp -o Exercise1_MPI
$ mpiexec -n 4 ./Exercise1_MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

double *generarArregloAleatorio(int numeroDeElementos);
double *generarArregloSen(int numeroDeElementos);

int main()
{
   int N = 10;
   int k;
   int i,iIni,iFin;
   double cpu_time_used;
   double *A,*B,*C;
   //==================
   // MPI [1]: definitions
   int numtasks, taskid;
   int Nx,Ndom;
   //==================

   //__________________________________________
   // Array of N elements
   A = generarArregloAleatorio(N);
   B = generarArregloAleatorio(N);
   C = (double *)malloc(N*sizeof(double));

   //==========================================
   // MPI [2]: start
   MPI_Init(NULL,NULL);                     /* Initialization */
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks); /* Number of processes */
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);   /* Identifier of each process */
   //==========================================

   iIni = 0;
   iFin = N;
   //==========================================
   // MPI [3]: Work division
   Ndom = (int)(1.0*N/numtasks);

   k = (int)(N % numtasks);
   if (taskid < k) {
       Ndom += 1;
       iIni = taskid*Ndom;
   }
   else {
       iIni = taskid*Ndom+k;
   }
   iFin = iIni + Ndom;
   //==========================================

   //__________________________________________
   // Operations: add two vectors
   clock_t start, end;
   start = clock();
   //-----------------------
   for (i = iIni; i < iFin; i++) {
        C[i] = A[i] + B[i];
        //printf("%f + %f = % f\n",A[i],B[i],C[i]);
   }
   //-----------------------
   end = clock();
   cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC;

   //__________________________________________
   // Print vector
   for (i=0; i<N; i++) {
        printf("%i : %f + %f = % f\n",i,A[i],B[i],C[i]);
    }


   //__________________________________________
   // Print
   printf("Processor %i of %i: \n \
           Ndom    = %i  \n \
           iIni  = %i  \n \
           iFin  = %i  \n \
           time  = %lf \n \n",taskid,numtasks,Ndom,iIni,iFin,cpu_time_used);

   //__________________________________________
   // Free memory
   free(A);
   free(B);
   free(C);

   //==========================================
   // MPI [4]: final
   MPI_Finalize();
   //==========================================

   return 0;
}

//sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
// Array of random numbers between 0 and 100

double *generarArregloAleatorio(int N){
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
// Integral of the sine

double *generarArregloSen(int N){
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
