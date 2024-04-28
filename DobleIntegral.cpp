/*sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!         PROGRAMA: Suma de los elementos de un vector          !
!                          PARALELO MPI                         !
!                    Miguel Angel Uh Zapata                     !
!                                                               !
!    HOW TO RUN :                                               !
!  $ mpicxx DobleIntegral.cpp -o run_DobleIntegral                  !
!  $ mpiexec -n 4 ./run_DobleIntegral                               !
!                                                               !
sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

double *generateRandomArray(int numberOfElements);
double *generateSinArray(int numberOfElements);

int main(int argc, char *argv[])
{
   int N = 100000;
   int M = 100000;
   int i,iIni,iFin;
   int j,jIni,jFin;
   double suma,cpu_time_used;
   double *x,*y; 
   double Dx,Dy;     

   //==================
   // MPI [1]: definiciones
   int divx = 2;
   int divy = 4;
   int numtasks, taskid;
   int Nx,Ndom,NdomF;
   int Ny,Mdom,MdomF;
   double sumaglob;
   //==================
   const int ndims = 2;
   int dims_vec[ndims]; 
   int periodicite[ndims]; 
   int reorganisation;
   int coords[ndims];
   MPI_Comm comm2D;    
   //==================
   
   //__________________________________________
   // Vector de N elementos
   
   x = (double *)malloc(N*sizeof(double));
   y = (double *)malloc(M*sizeof(double));
   
   Dx = 2.0/(N-1);   
   for (int i = 0; i < N; i++){
        x[i] = 0.0 + i*Dx;
        //printf("%d %f\n",i,x[i]);
   }

   Dy = 2.0/(M-1);   
   for (int j = 0; j < M; j++){
        y[j] = -1.0 + j*Dy;
        //printf("%d %f\n",j,y[j]);
   }
            
   //==========================================
   // MPI [2]: inicio
   MPI_Init(NULL,NULL);                     /* Initializacion*/
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks); /* Numero de procesos */
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);   /* Identificador de cada proceso */
   //==========================================

   if (divx*divy==numtasks){
      printf("\n CORRECTO: numtasks = divx*divy = %d \n",divx*divy);
      }
   else {
      printf("\n INCORRECTO: numtasks = %d  y  divx*divy = %d \n",
      numtasks,divx*divy);
      return 0;
   } 
 
   //==========================================
   // MPI [2]: Topologia (cartesiana)
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
   printf("\ntaskid=%d coords (%d,%d) \n",taskid,coords[0],coords[1]);
      
   //==========================================
   // MPI [3]: Division del trabajo
   iIni = 0;
   iFin = N;   
   //-----------------
   Ndom  = (int)(1.0*N/divx);
   NdomF = N - Ndom*(divx-1);
   //-----------------
   Nx = Ndom ;
   if (coords[0]==divx-1){Nx = NdomF-1;}
   //-----------------
   iIni = coords[0]*Ndom;
   iFin = iIni + Nx;
   //==========================================

   //==========================================
   // MPI [3]: Division del trabajo
   jIni = 0;
   jFin = M;   
   //-----------------
   Mdom  = (int)(1.0*M/divy);
   MdomF = M - Mdom*(divy-1);
   //-----------------
   Ny = Mdom ;
   if (coords[1]==divy-1){Ny = MdomF-1;}
   //-----------------
   jIni = coords[1]*Mdom;
   jFin = jIni + Ny;
   //==========================================
       
   //__________________________________________
   // Operaciones: sumar sus componentes
   clock_t start, end;
   start = clock();
   
   //----------------------- 
   suma = 0.0;
   for (i = iIni; i < iFin; i++) {
       for (j = jIni; j < jFin; j++) {
           suma = suma + (1.0 + 6.0*x[i]*y[j]*y[j])*Dx*Dy;
       } 
   }
   //-----------------------

   //==========================================
   // MPI [4]: sumar componentes de cada procesador
   MPI_Allreduce(&suma,&sumaglob,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2D);
   //==========================================

   end = clock();
   cpu_time_used = ((double)(end-start))/CLOCKS_PER_SEC;
   
   //__________________________________________
   // Imprimir
   printf("Procesor %i de %i: \n \
           Nx    = %i  \n \
           iIni  = %i  \n \
           iFin  = %i  \n \
           Ny    = %i  \n \
           jIni  = %i  \n \
           jFin  = %i  \n \
           suma  = %f  \n \
           sumaT = %f  \n \
           Error = %f  \n \
           time  = %lf \n \n",taskid,numtasks,
           Nx,iIni,iFin,
           Ny,jIni,jFin,
           suma,sumaglob,
           12.0 - sumaglob,
           cpu_time_used);

   //__________________________________________
   // Liberar memoria
   free(x);
   free(y);

   //==========================================
   // MPI [5]: final
   MPI_Finalize();
   //==========================================

   return 0;
}



//sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss