!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!---------------------------------------------------------------------!
!          PROGRAMA: Solucion de la ecuacion de calor con MPI         !
!                       Miguel Angel Uh Zapata                        !
!                               2024                                  !
!---------------------------------------------------------------------!
!     HOW TO RUN:                                                     !
!     $ mpif77 Calor2DMPI2024_Opt1.F90 -o run_Calor2DMPI                  !
!     $ mpiexec -n 4 ./Calor2DMPI                                     !
!---------------------------------------------------------------------!
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      program Calor2DMPI2024_Opt1
      
      implicit none
      include 'mpif.h'       !<--- MPI paqueteria

      integer :: i,j,n,index
      integer :: Nx,Ny,Nt
      integer :: NxG,NyG,NN,NNF,nrec,iIni,iFin
      real*8, parameter :: pi = 3.14159265359d0 
      real*8  :: kx,ky,xI,xF,yI,yF,tF,tI,Dx,Dy,Dt,t
      real*8  :: rx,ry,aC,aN,aS,aE,aW
      real*8  :: suma,sumaglob
      real    :: inicio,final,tiempo
!     -------------------
      real*8, dimension(:),  allocatable :: x,y
      real*8, dimension(:),  allocatable :: xG,yG
      real*8, dimension(:,:),allocatable :: Uold,Unew
      real*8, dimension(:,:),allocatable :: UG,PG
      
!     .----------------------------------------------------.
!     |                  MPI: definiciones                 |
!     .----------------------------------------------------.

      integer :: ierr
      integer :: numtasks,taskid
      integer :: tipo_col,tipo_block
      integer, dimension(2) :: vecino
      integer, dimension(MPI_STATUS_SIZE) :: statut       
      integer, dimension(:),  allocatable :: index_global            
      integer :: etiqueta
      etiqueta = 2001
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                           INICIALIZACION                            !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

!     ______________________________________________________
!     Parametros 
 
      kx = 1.0d0        ! Termino difusivo en x
      ky = 1.0d0        ! Termino difusivo en y
!     ----------
      xI = 0.0d0        ! Inicio del dominio en x
      xF = 1.0d0        ! Fin del dominio en x
      yI = 0.0d0        ! Inicio del dominio en y
      yF = 1.0d0        ! Fin del dominio en y
!     ----------
      tI = 0.0d0        ! Tiempo inicial
      tF = 0.2d0        ! Tiempo final
!     ----------      
      Nt  = 200000      ! Numero de puntos en t
      NxG = 150         ! Numero de puntos en x (GLOBAL)
      NyG = 150         ! Numero de puntos en y (GLOBAL)

!     ______________________________________________________
!     Discretización

      Dx = (xF-xI)/(NxG-1) 
      Dy = (yF-yI)/(NyG-1)  
      Dt = (tF-tI)/(Nt-1)      
      rx  = kx*(Dt/(Dx**2))
      ry  = ky*(Dt/(Dy**2))

!     =========================================          
!     Variables globales

      allocate(xG(NxG))
      allocate(yG(NyG))
      allocate(UG(NxG,NyG))
      allocate(PG(NxG,NyG))

!     ------------------------
!     Malla
      do i=1,NxG
         xG(i) = xI + (i-1)*Dx
      enddo
      do j=1,NyG
         yG(j) = yI + (j-1)*Dy
      enddo

!     ------------------------
!     Condiciones Iniciales
      do i=1,NxG
         do j=1,NyG
            UG(i,j) = 3.0d0*sin(pi*xG(i)+pi*yG(j))**2
            !UG(i,j) = sin(xG(i)+yG(j))**2
         enddo
      enddo

!     ------------------------
!     Condiciones de frontera
      do j=1,NyG
         UG(1,j)   = 2.0d0  !Oeste
         UG(NxG,j) = 1.0d0  !Este
      enddo
      do i=1,NxG
         UG(i,1)   = 1.0d0  !Sur
         UG(i,NyG) = 3.0d0  !Norte
      enddo

!     ------------------------
      PG = 0.0d0
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                         COMPONENTES MPI                             !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
      
!     .----------------------------------------------------.
!     |                  MPI: Initialization               |
!     .----------------------------------------------------.
   
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)     
      call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)

!     .----------------------------------------------------.
!     |             MPI: Topologia (cartesiana)            |
!     .----------------------------------------------------.

      vecino(1) = MPI_PROC_NULL
      vecino(2) = MPI_PROC_NULL
      if (taskid.gt.0)          vecino(1) = taskid - 1
      if (taskid.lt.numtasks-1) vecino(2) = taskid + 1
      
!     .----------------------------------------------------.
!     |             MPI: División del dominio              |
!     .----------------------------------------------------.

!     =========================================
!     División del dominio en x
      NN  = floor(1.0d0*NxG/numtasks)
      NNF = NxG - NN*taskid
!     -------------      
      Nx  = NN + 2
      if (taskid.eq.0)          Nx = NN  + 1
      if (taskid.eq.numtasks-1) Nx = NNF + 1
!     =========================================
!     División del dominio en y      
      Ny = NyG
      
!     .----------------------------------------------------.
!     |           MPI: Indices locales y globales          |
!     .----------------------------------------------------.

      allocate(index_global(Nx))
      
      if (taskid.eq.0) then
         do i=1,Nx
            index_global(i) = i
         enddo
      else
         do i=1,Nx
            index_global(i) = taskid*NN + i-1 
         enddo
      endif

!     -------------
!     IMPRIMIR
!     ------------- 
      if (taskid.eq.0) then
         print*,''
         print*,'(NxG,NyG)=',NxG,NyG
         print*,'Criterio CFL < 1/2 :',rx+ry
         print*,''
      endif
                  
      write(*,'(a9,i3,a6,i4,a12,i4,a16,i8,i8)') &
             'taskid =',taskid, &
             ', NN =',NN, &
             ', Nx_local =',Nx, &
             ', index global =', index_global(1),index_global(Nx)      

!     .----------------------------------------------------.
!     |               MPI: Variables locales               |
!     .----------------------------------------------------.
      
      allocate(x(Nx))
      allocate(y(Ny))
      allocate(Uold(Nx,Ny))
      allocate(Unew(Nx,Ny))

!     ------------------------
!     Malla
      do i=1,Nx
         index = index_global(i)
         x(i) = xG(index)
      enddo
      do j=1,Ny
         y(j) = yG(j)
      enddo

!     ------------------------
!     Condiciones Iniciales
      do i=1,Nx
         index = index_global(i)
         do j=1,Ny
            Uold(i,j) = UG(index,j)
         enddo
      enddo

!     ------------------------
!     Actualizacion
      Unew = Uold
             
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                        PROGRAMA PRINCIPAL                           !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

!     ______________________________________________________          
!     Guardar condiciones iniciales y de frontera
      nrec = 0
      if (taskid.eq.0) call save_paraview(NxG,NyG,xG,yG,t,UG,PG,nrec)
      
!     .----------------------------------------------------.
!     |        MPI: Tipo de vectores (comunicaciones)      |
!     .----------------------------------------------------.
!           vector of Ny elements separated by Nx

      call MPI_TYPE_VECTOR(Ny,1,Nx,MPI_DOUBLE_PRECISION,&
                           tipo_col,ierr)
      call MPI_TYPE_VECTOR(1,NN*Ny,1,MPI_DOUBLE_PRECISION,&
                           tipo_block,ierr)
      call MPI_TYPE_COMMIT(tipo_col,ierr)
      call MPI_TYPE_COMMIT(tipo_block,ierr)
      

!     .----------------------------------------------------.
!     |               MPI: Loop de calculos                |
!     .----------------------------------------------------.

      call cpu_time(inicio)
      
      aC = 1-2.0d0*(rx+ry)
      aE = rx
      aW = rx
      aS = ry
      aN = ry

      DO n=1,Nt
!        ----------------------------------------------------
!        Nuevos valores
         do i=2,Nx-1
            do j=2,Ny-1
               Unew(i,j) = aC*Uold(i,j)   &
                         + aW*Uold(i-1,j) &
                         + aE*Uold(i+1,j) &
                         + aS*Uold(i,j-1) &
                         + aN*Uold(i,j+1) 
            enddo
         enddo
         
!        ----------------------------------------------------
!        Comunicacion
         call MPI_SENDRECV(Unew(2,1) ,1,tipo_col,vecino(1),etiqueta, &
                           Unew(Nx,1),1,tipo_col,vecino(2),etiqueta, &
                           MPI_COMM_WORLD,statut,ierr)

         call MPI_SENDRECV(Unew(Nx-1,1),1,tipo_col,vecino(2),etiqueta,&
                           Unew(1,1)   ,1,tipo_col,vecino(1),etiqueta,&
                           MPI_COMM_WORLD,statut,ierr)

!        ----------------------------------------------------
!        Actualización
         Uold = Unew
         
!        ----------------------------------------------------
!        Guardar resultados
         if (mod(n,1000).eq.0) then   
!           ----------------------- 
            t = tI + n*Dt
            nrec = nrec + 1
            if (taskid.eq.numtasks-1) print*,nrec,' time=',t 
!           -----------------------
            IF (NxG.le.150) THEN
            if (taskid.lt.numtasks-1) then
               iIni = 2
               if (taskid.eq.0) iIni = 1
               iFin = iIni + NN-1
               call MPI_SEND(Unew(Iini:iFin,1:Ny),1,tipo_block, &
                             numtasks-1,etiqueta,MPI_COMM_WORLD,ierr)
            endif  
            
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
                            
            if (taskid.eq.numtasks-1) then
               iIni = taskid*NN + 1
               iFin = iIni + Nx
               UG(iIni:iFin,1:Ny) = Unew(2:Nx,1:Ny)  
               PG = numtasks-1    
               do i=0,numtasks-2
                  iIni = i*NN + 1
                  iFin = iIni + NN-1
                  call MPI_RECV(UG(iIni:IFin,1:Ny),1,tipo_block, &
                                i,etiqueta,MPI_COMM_WORLD,statut,ierr)
                  PG(iIni:iFin,1:Ny) = i              
               enddo
               call save_paraview(NxG,NyG,xG,yG,t,UG,PG,nrec)
            endif
!           -----------------------
            ENDIF         
         endif 
           
      ENDDO

      call cpu_time(final)

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                           VERIFICAR RESULTADOS                      !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

!     ______________________________________________________
!     Imprimir suma de los elementos para verificar
      
      iIni = 2
      iFin = Nx-1 
      if (taskid.eq.0)          iIni = 1
      if (taskid.eq.numtasks-1) iFin = Nx 
!     -----------------------       
      suma = 0.0d0
      do i=iIni,iFin
         do j=1,Ny
            suma = suma + Unew(i,j)
         enddo
      enddo
!     -----------------------
      call MPI_ALLREDUCE(suma,sumaglob,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,MPI_COMM_WORLD,ierr)
      suma = sumaglob
!     -----------------------
      print*,'taskid=',taskid,', Suma total = ',suma

!     ______________________________________________________
!     Tiempo final
      tiempo = final-inicio   
      print*,'tiempo=',tiempo,', taskid=',taskid     
 
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                             FINALIZACION                            !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

!     .----------------------------------------------------.
!     |                    MPI: Completion                 |
!     .----------------------------------------------------.

      call MPI_TYPE_FREE(tipo_col,ierr)
      call MPI_TYPE_FREE(tipo_block,ierr)
      call MPI_FINALIZE(ierr) 

!     _________________________________________    
!     Liberar
      deallocate(x,y,Uold,Unew)      
      deallocate(xG,yG,UG,PG,index_global)
    
      return
      end

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                   GUARDANDO RESULTADOS EN PARAVIEW                  !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      SUBROUTINE save_paraview(Nx,Ny,x,y,t,Unew,PG,nrec)
      
      implicit none     
      
      integer :: Nx,Ny,i,j,nrec,irec
      real*8  :: t
      real*8, dimension(:)   :: x(Nx),y(Ny)
      real*8, dimension(:,:) :: Unew(Nx,Ny),PG(Nx,Ny)
      character*27 filen
      
!      ________________________________________________________
!     |                                                        |
!     |                         Formats                        |
!     |________________________________________________________|

 2    format(1(1x,e12.5))
 3    format(3(1x,e12.5))

 100  format('# vtk DataFile Version 2.0')
 110  format('Sample rectilinear grid')
 120  format('ASCII')
 130  format('DATASET RECTILINEAR_GRID')
 140  format('DIMENSIONS ',i5,1x,i5,' 1')
 150  format('X_COORDINATES ',i5,' float')
 160  format('Y_COORDINATES ',i5,' float')
 170  format('Z_COORDINATES 1 float')
 180  format('0.0')
 190  format('POINT_DATA ',i7)

!      ________________________________________________________
!     |                                                        |
!     |                       File name                        |
!     |________________________________________________________|

      irec=60
      filen='P-    .vtk'
      write(filen(3:6),'(i4.4)') nrec

!      ________________________________________________________
!     |                                                        |
!     |                          Save                          |
!     |________________________________________________________|
      
!     ___________________________________________________
!     Open
      open(irec,FILE=filen)
!     ___________________________________________________
!     Title head
      write(irec,100)
      write(irec,110)
      write(irec,120)
      write(irec,130)
      write(irec,140) Nx,Ny
!     ___________________________________________________
!     Coordinates
      write(irec,150) Nx
      do i=1,Nx
         write(irec,2) x(i)
      enddo
      write(irec,160) Ny
      do j=1,ny
         write(irec,2) y(j)
      enddo
      write(irec,170)
      write(irec,180)
!     ___________________________________________________
!     Number of points
      write(irec,190) Nx*Ny
!     ___________________________________________________
!     Scalars
      write(irec,*) 'SCALARS U float'
      write(irec,*) 'LOOKUP_TABLE default'
      do j=1,Ny
         do i=1,Nx
            write(irec,2) Unew(i,j) 
         enddo
      enddo
      write(irec,*) 'SCALARS Proc float'
      write(irec,*) 'LOOKUP_TABLE default'
      do j=1,Ny
         do i=1,Nx
            write(irec,2) PG(i,j) 
         enddo
      enddo    
!     ___________________________________________________
!     Close
      rewind(irec)
      close(irec)
         
!      ________________________________________________________
!     |                                                        |
!     |                       Display                          |
!     |________________________________________________________|

      write(*,*) '  '
      write(*,'(t5,a12,i4,a9,g10.3)') 'Paraview No.',nrec,', time :',t       

      return
      end
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                         FIN DEL PROGRAMA                            !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!