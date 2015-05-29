!Ising Model
program MAIN
implicit none
	integer, allocatable, DIMENSION(:,:) :: arr
	real :: num,Temp,timestart,timeend,time,Eavg,Mavg
	integer*4 N, i,j,t,E,E0,M,E1,E2
	Eavg = 0
	t= 100
	open (unit =11, file = "IsingEavg.txt") 
	open (unit =12, file = "IsingMavg.txt")
	open (unit =13, file = "timef90.txt")



	do N = 10, 52, 2
		ALLOCATE(arr(N,N))
		call cpu_time(timestart)
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
		do j = 1, 51, 1
			temp = j/10.
			Eavg= 0
			Mavg = 0
			!WARM UP
			call flip(arr,E0,E1,N,temp)
			call calcE(N,arr,E1,M)
			call flip(arr,E0,E2,N,temp)
			call calcE(N,arr,E2,M)
			do while (E1-E2>N)
				!Choose a an Energy 3 iterations down and checks its 					difference against initial energy	
				E1=E2 
				call flip(arr,E0,E2,N,temp)
				call calcE(N,arr,E2,M)
				call flip(arr,E0,E2,N,temp)
				call calcE(N,arr,E2,M)
				call flip(arr,E0,E2,N,temp)
				call calcE(N,arr,E2,M)
			enddo
			!End Warm Up
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
		write(13,*), time
		DEALLOCATE(arr)
	end do
	print *, time

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
