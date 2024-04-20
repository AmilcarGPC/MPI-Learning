program main
      ! Declaration of variables & functions
      include 'mpif.h'
      double precision  PI25DT
      parameter        (PI25DT = 3.141592653589793238462643d0)
      double precision  mypi, pi, h, sum, x, f, a
      integer n, myid, numprocs, i, ierr
      real start, finish
      double precision :: min_time, max_time, total_time
  
      ! Function to integrate
      f(a) = 4.d0 / (1.d0 + a*a)
  
      ! MPI: Initialization
      call MPI_INIT(ierr)
      ! MPI: Assign ID to each processor
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      ! MPI: Read number of processors
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  
      ! Processor 0: read n
      if (myid.eq.0) then
          print*, 'Enter the number of intervals:  '
          read(*,*) n
      endif
  
      ! MPI: Broadcast n
      call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  
      call cpu_time(start)
  
      ! Calculate the interval size
      h = 1.0d0 / n
      sum  = 0.0d0
      do i = myid + 1, n, numprocs
          x = h * (dble(i) - 0.5d0)
          sum = sum + f(x)
      enddo
      mypi = h * sum
  
      ! MPI: Collect all the partial sums
      call MPI_REDUCE(mypi, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
      call cpu_time(finish)
  
      ! Processor 0: print the answer and error
      if (myid.eq.0) then
          print*, 'pi = ', pi, ' Error =', abs(pi - PI25DT)
      endif
  
      ! Calculate CPU time for each process
      cpu_time_used = finish - start
  
      ! MPI: Reduction for CPU times
      call MPI_REDUCE(cpu_time_used, min_time, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(cpu_time_used, max_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(cpu_time_used, total_time, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
      ! Processor 0: print minimum, maximum, and average CPU times
      if (myid.eq.0) then
          print*, 'Min CPU time =', min_time, 'sec'
          print*, 'Max CPU time =', max_time, 'sec'
          print*, 'Average CPU time =', total_time / numprocs, 'sec'
      endif
  
      ! MPI: Finalization
      call MPI_FINALIZE(ierr)
  
  end
  