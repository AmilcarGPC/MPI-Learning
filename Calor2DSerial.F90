!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!---------------------------------------------------------------------!
!          PROGRAMA: Solucion de la ecuacion de calor en Serie        !
!                       Miguel Angel Uh Zapata                        !
!                               2024                                  !
!---------------------------------------------------------------------!
!    HOW TO RUN:                                                      !
!    $ mpif77 Calor2DSerial.F90 -o run_Calor2DSerial                      !
!    $ ./run_Calor2DSerial                                                !
!---------------------------------------------------------------------!
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      program Calor2DSerial
      
      implicit none     
      
      integer :: i,j,n
      integer :: Nx,Ny,Nt
      real*8, parameter :: pi = 3.14159265359d0 
      real*8  :: kx,ky      
      real*8  :: xI,xF,yI,yF
      real*8  :: tF,tI
      real*8  :: Dx,Dy,Dt,t
!     -------------------
      real*8, dimension(:),  allocatable :: x
      real*8, dimension(:),  allocatable :: y
      real*8, dimension(:,:),allocatable :: Uold
      real*8, dimension(:,:),allocatable :: Unew
!     -------------------
      real*8  :: rx,ry
      real*8  :: aC,aN,aS,aE,aW
!     -------------------
      real*8  :: suma
      real    :: inicio,final,tiempo
!     -------------------
      integer :: nrec


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
      Nt = 200000       ! Numero de divisiones en el tiempo
      Nx = 400          ! Numero de puntos en x
      Ny = 400          ! Numero de puntos en y

!     ______________________________________________________
!     Discretización

      Dx = (xF-xI)/(Nx-1) 
      Dy = (yF-yI)/(Ny-1)  
      Dt = (tF-tI)/(Nt-1)      
      rx  = kx*(Dt/(Dx**2))
      ry  = ky*(Dt/(Dy**2))
!     -------------
!     IMPRIMIR      
      print*,'(Nx,Ny)=',Nx,Ny
      print*,'Criterio CFL < 1/2 :',rx+ry
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                               MALLAS                                !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      allocate(x(Nx))
      allocate(y(Ny))

      do i=1,Nx
         x(i) = xI + (i-1)*Dx
      enddo
      do j=1,Ny
         y(j) = yI + (j-1)*Dy
      enddo
       
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                        PROGRAMA PRINCIPAL                           !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      call cpu_time(inicio)
      
!     ______________________________________________________          
!     Condiciones en la solución

      allocate(Uold(Nx,Ny))
      allocate(Unew(Nx,Ny))
      
!     ------------------------
!     Condiciones iniciales
      do i=1,Nx
         do j=1,Ny
            Uold(i,j) = 3.0d0*sin(pi*x(i)+pi*y(j))**2
         enddo
      enddo
      Unew = Uold

!     ------------------------
!     Condiciones de frontera
      do j=1,Ny
         Unew(1,j)  = 2.0d0  !Oeste
         Unew(Nx,j) = 1.0d0  !Este
      enddo
      do i=1,Nx
         Unew(i,1)  = 1.0d0  !Sur
         Unew(i,Ny) = 3.0d0  !Norte
      enddo

      nrec = 0
      call save_paraview(Nx,Ny,x,y,t,Unew,nrec)
      
!     ______________________________________________________          
!     Calculos en los puntos interiores

      aC = 1-2.0d0*(rx+ry)
      aE = rx
      aW = rx
      aS = ry
      aN = ry
               
      DO n=1,Nt
!        ------------------------
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
!        ------------------------
!        Actualización
         Uold = Unew
      ENDDO

      call cpu_time(final)
      tiempo = final-inicio
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                           VERIFICAR RESULTADOS                      !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

!     ______________________________________________________
!     Imprimir suma de los elementos para verificar
            
      suma = 0.0d0
      do i=1,Nx
         do j=1,Ny
            suma = suma + Unew(i,j)
         enddo
      enddo
!     -----------------
!     IMPRIMIR: suma 
      print*, '  ' 
      print*,'Suma total = ',suma
        
!     ______________________________________________________
!     IMPRIMIR: tiempo final 
      print*,'tiempo=',tiempo
      print*, '  '

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                             FINALIZACION                            !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      deallocate(x)
      deallocate(y)
      deallocate(Uold)
      deallocate(Unew)            
    
      return
      end

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!                                                                     !
!                   GUARDANDO RESULTADOS EN PARAVIEW                  !
!                                                                     !
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      SUBROUTINE save_paraview(Nx,Ny,x,y,t,Unew,nrec)
      
      implicit none     
      
      integer :: Nx,Ny,i,j,nrec,irec
      real*8  :: t
      real*8, dimension(:)   :: x(Nx),y(Ny)
      real*8, dimension(:,:) :: Unew(Nx,Ny)
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
      filen='U-    .vtk'
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