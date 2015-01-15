!Ising Model
program MAIN
implicit none
	integer, allocatable, DIMENSION(:,:) :: arr
	real :: num,Temp,timestart,timeend,time,Eavg,Mavg
	integer*4 N, i,j,t,E,E0,M
	Eavg = 0
	open (unit =11, file = "IsingEavg.txt") 
	open (unit =12, file = "IsingMavg.txt")

	call cpu_time(timestart)

	write(*,*) 'Enter the size of matrix: '
	read(*,*) N
	write(*,*) 'Enter the amount of iterations: '
	read(*,*) t
	ALLOCATE(arr(N,N))


	do j = 1, N
		do i = 1,N
			CALL RANDOM_NUMBER(num)
			if(num<0.5) Then
				arr(i,j) = 1
			Else
				arr(i,j) = -1
			END IF
		end do
	end do
	print *, arr
	do j = 1, 51, 1
		temp = j/10.
		Eavg= 0
		Mavg = 0
		do i = 1, t
			call flip(arr,E0,E,N,temp)
			call calcE(N,arr,E,M)
			Eavg = E+Eavg
			Mavg = M+Mavg

		end do
		Eavg = Eavg/(N*N*t)
		Mavg = Mavg/(N*N*t)
		write(11,*), Eavg
		write(12,*), Mavg
	end do
	call cpu_time(timeend)
        time = timeend-timestart
	print *, time

	DEALLOCATE(arr)
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
		arrnew = arr
		call localE(i,j,N,arr,Enew)
		E0 = Enew
		if(arr(i,j) == 1) Then
			arrnew(i,j) = -1
			call localE(i,j,N,arrnew,Enew)
		Else
			arrnew(i,j) = 1
			call localE(i,j,N,arrnew,Enew)
		END IF

		if(Enew<=E0) Then
			arr = arrnew
		else
			p = exp(-(Enew-E0)/temp)
			CALL RANDOM_NUMBER(num)
			if(p>=num) Then
				arr = arrnew
			END IF
		END IF

	60 continue
  50  continue

	end subroutine flip

!***********************************************************************
