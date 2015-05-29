!Ising Model
!To compile fortran code with MPI
!mpif90 Isingmpi.f90 -o Isingmpi
!Run the MPI code with 4 proccessors
!mpirun -np 4 ./Isingmpi

program MAIN

  	! Include the MPI library definitons:
 	 include 'mpif.h'
 	integer numtasks, rank, ierr, rc, len, i
  	character*(MPI_MAX_PROCESSOR_NAME) name

	integer, allocatable, DIMENSION(:,:) :: arr
	real, DIMENSION(50) :: Eavg,Mavg
	real :: num,Temp,timestart,timeend,time
	integer*4 N, k,j,t,E,E0,M,x
	Eavg = 0
	Mavg = 0
	t= 100


	  ! Initialize the MPI library:
  call MPI_INIT(ierr)
  if (ierr .ne. MPI_SUCCESS) then
     print *,'Error starting MPI program. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  end if

  ! Get the number of processors this job is using:
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

  ! Get the rank of the processor this thread is running on.  (Each
  ! processor has a unique rank.)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  ! Get the name of this processor (usually the hostname)
  call MPI_GET_PROCESSOR_NAME(name, len, ierr)
  if (ierr .ne. MPI_SUCCESS) then
     print *,'Error getting processor name. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  end if
	open (unit =11, file = "IsingEavg1.txt")
	open (unit =12, file = "IsingEavg2.txt") 
	open (unit =13, file = "IsingEavg3.txt") 
	open (unit =14, file = "IsingEavg4.txt") 
	open (unit =15, file = "IsingMavg1.txt")
	open (unit =16, file = "IsingMavg2.txt")
	open (unit =17, file = "IsingMavg3.txt")
	open (unit =18, file = "IsingMavg4.txt")
	open (unit =19, file = "timef90mpi.txt")

	call cpu_time(timestart)

	!call MPI_Barrier(MPI_COMM_WORLD)
	do N = 10, 52, 2
		call cpu_time(timestart)
		ALLOCATE(arr(N,N))

		if(rank==0)then
			do j = 1, N
				do k = 1,N
					CALL RANDOM_NUMBER(num)
					if(num<0.5) Then
						arr(k,j) = 1
					Else
						arr(k,j) = -1
					END IF
				end do
			end do

		end if



		temp = rank/10.

		do while(temp<5.0)
			x = int(temp*10)
			Eavg(x)= 0
			Mavg(x) = 0
			do k = 1, t
				call flip(arr,E0,E,N,temp)
				call calcE(N,arr,E,M)
				Eavg(x) = E+Eavg(x)
				Mavg(x) = M+Mavg(x)

			end do
			Eavg(x) = Eavg(x)/(N*N*t)
			Mavg(x) = Mavg(x)/(N*N*t)
			temp = temp + numtasks/10.

		end do


	
		write(11+rank,*), Eavg
		write(15+rank,*), Mavg
		
		call cpu_time(timeend)
        	time = timeend-timestart
		if(rank==2)then
			write(19,*), time
		end if
		DEALLOCATE(arr)

	end do
	! Tell the MPI library to release all resources it is using:
  	call MPI_FINALIZE(ierr)
end program Main

!***********************************************************************
      subroutine calcE(N,arr,E,M)
	integer, INTENT(INOUT), DIMENSION(N,N) :: arr
	integer, INTENT(INOUT) :: E,M
	integer :: i,j
	E = 0
	M = 0
        do 30 i = 1, N
         do 40 j = 1, N
		M = M + arr(i,j)
		if(i<N) Then
			E=E -arr(i,j)*arr(i+1,j)
		else if(i==N) Then
			E=E -arr(i,j)*arr(1,j)
		END IF

		if(j<N) Then
			E=E -arr(i,j)*arr(i,j+1)
		else if(j==N) Then
			E=E -arr(i,j)*arr(i,1)
		END IF
	 40 continue
  	30 continue
	end subroutine calcE

!***********************************************************************
!Calculates the local energy 
      subroutine localE(i,j,N,arr, Enew)
	integer, INTENT(INOUT), DIMENSION(N,N) :: arr
	integer, INTENT(INOUT) :: Enew
	integer, INTENT(IN) :: i,j
	E = 0
	if(i>1 .and. i<N) Then
		E=E -arr(i,j)*arr(i+1,j)-arr(i,j)*arr(i-1,j)
	else if(i==N) Then
		E=E -arr(i,j)*arr(1,j)-arr(i,j)*arr(i-1,j)
	else if(i==1) Then
		E=E -arr(i,j)*arr(i+1,j)-arr(i,j)*arr(N,j)
	END IF
	if(j>1 .and. j<N) Then
		E=E -arr(i,j)*arr(i,j+1)-arr(i,j)*arr(i,j-1)
	else if(j==N) Then
		E=E -arr(i,j)*arr(i,1)-arr(i,j)*arr(i,j-1)
	else if(j==1) Then
		E=E -arr(i,j)*arr(i,j+1)-arr(i,j)*arr(i,N)
	END IF
	Enew = E
	end subroutine localE

!***********************************************************************
 !Takes the array and sweeps through it flipping each value, and checking the old energy vs the new energy     

	subroutine flip(arr,E0,E,N,temp)
	integer, INTENT(INOUT), DIMENSION(N,N) :: arr
	integer, DIMENSION(N,N) :: arrnew
	real :: p,num, temp 
	integer :: i, j, Enew,E,E0
	do 50 j = 1, N
         do 60 i = 1, N
		call localE(i,j,N,arr,Enew)
		E0 = Enew
		Enew = Enew*(-1)
		if(Enew<=E0) Then
			arr(i,j) = (-1)*arr(i,j)
		else
			p = exp(-(Enew-E0)/temp)
			CALL RANDOM_NUMBER(num)
			if(p>=num) Then
				arr(i,j) = (-1)*arr(i,j)
			END IF
		END IF

	60 continue
  50  continue

	end subroutine flip

!***********************************************************************
