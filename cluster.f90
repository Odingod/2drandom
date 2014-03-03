PROGRAM bd
  ! These are comments
  !
  ! Ballistic Deposition (BD) and RD simulations compared.
  !
  ! Compile as: gfortran -o BD_and_RD BD_and_RD.f90
  !
  IMPLICIT NONE ! This forces us to define all variables
  INTEGER, ALLOCATABLE ::  h(:) ! height
  REAL(KIND(1.d0)), ALLOCATABLE ::  w(:) ! roughness(^2)
  INTEGER :: i, j, N, Nh, Nr, ii, jj, kk
  REAL :: ra
  CHARACTER(len=40) :: w_file
  LOGICAL :: do_bd=.TRUE.

  CALL init_random_seed()

  ! WRITE(*,*) 'Do cluster growth? (false makes RD), say t or f '
  ! READ(*,*) do_bd
  ! WRITE(*,*) 'Number of sites (horizontal size of cluster) [try 100]'
  ! READ(*,*) N
  ! ALLOCATE(h(0:N+1))
  ! WRITE(*,*) 'How many particles (per column)? [try 1000]'
  ! READ(*,*) Nh
  ! WRITE(*,*) 'How many runs? [try 100]'
  ! READ(*,*) Nr
  ! ALLOCATE(w(Nh))

  w=0.d0
  Nr = 200
  Nh = 10000
  DO kk=4, 9 
	  N = 2**kk
	  IF (ALLOCATED(h)) DEALLOCATE(h)
	  IF (ALLOCATED(w)) DEALLOCATE(w)
	  ALLOCATE(w(Nh))
	  ALLOCATE(h(0:N+1))
	  WRITE(*,*) 'N on ', N
	  w=0.d0
	  IF (do_bd) THEN 
		  WRITE(w_file,'(a,i0,a)') 'w_average', N, '.txt'
		  ELSE
			 WRITE(w_file,'(a,i0,a)') 'rd_w_average', N, '.txt'     
		  END IF
	  OPEN(123,file=w_file)
	  DO jj=1, Nr
		 WRITE(*,*) 'Run ', jj, '/', Nr, kk-3, '/6' 
		 h=0  ! This contains heights.
		 DO ii=1, Nh
			DO j=1, N
			   i=random_site()
			   IF (do_bd) THEN
					IF (h(i-1) < h(i) .AND. h(i+1) < h(i)) THEN
						CALL RANDOM_NUMBER(ra)
						IF (ra  < 0.5) THEN
							h(i-1) = h(i-1)+1
							IF (i==1) h(N)=h(0)
							IF (i==2) h(N+1)=h(1)
						ELSE 
							h(i+1) = h(i+1)+1
							IF (i==N) h(1)=h(N+1)
							IF (i==N-1) h(0)=h(N)
						END IF
					ELSE IF (h(i-1)<h(i)) THEN
						h(i-1) = h(i-1) + 1
						IF (i==1) h(N)=h(0)
						IF (i==2) h(N+1)=h(1)
					ELSE IF (h(i+1)<h(i)) THEN
						h(i+1) = h(i+1) + 1
						IF (i==N) h(1)=h(N+1)
						IF (i==N-1) h(0)=h(N)
					ELSE 
						h(i) = h(i) + 1
						IF (i==1) h(N+1)=h(1) ! These take care of 
						IF (i==N) h(0)=h(N)
					END IF
					 ! the periodic boundary conditions.
			   ELSE
				  h(i)=h(i)+1 ! This is the rule of RD.
			   END IF
			END DO
			w(ii)=w(ii)+rough2(h(1:N))
		 END DO
	  END DO
	  w=w/REAL(Nr,KIND(1.d0))
	  DO ii=1, Nh
		 WRITE(123,*) ii, SQRT(w(ii))
	  END DO
	  CLOSE(123)
	  WRITE(*,*) 'Data in ', w_file
	END DO
CONTAINS
  FUNCTION random_site()
	INTEGER :: random_site
	REAL(KIND(1.d0)) :: r
	CALL RANDOM_NUMBER(r)
	random_site=MIN(INT(r*N)+1,N)
  END FUNCTION random_site
  FUNCTION rough2(h)
	REAL(KIND(1.d0)) :: rough2, eh, eh2
	INTEGER :: h(:), N
	REAL(KIND(1.d0)), DIMENSION(SIZE(h)) :: nh
	N=SIZE(h)
	nh=h-MINVAL(h)
	eh=SUM(nh)/REAL(N,KIND(1.d0))
	eh2=SUM(nh**2)/REAL(N,KIND(1.d0))
	rough2=eh2-eh**2
  END FUNCTION rough2
  SUBROUTINE init_random_seed()
	! This is the initialization of the random number generator
	IMPLICIT NONE
	INTEGER, ALLOCATABLE :: seed(:)
	INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
	INTEGER(8) :: count, tms

	CALL RANDOM_SEED(size = n)
	ALLOCATE(seed(n))
	WRITE(*,'(a,i0,a)') 'Random number generator seed has ', n,' numbers.'
	! First try if the OS provides a random number generator
	OPEN(newunit=un, file="/dev/urandom", access="stream", &
		 form="unformatted", action="read", status="old", iostat=istat)
	IF (istat == 0) THEN
	   WRITE(*,*) 'Found a stream, good.'
	   READ(un) seed
	   CLOSE(un)
	ELSE
	   ! Fallback to XOR:ing the current time and pid. The PID is
	   ! useful in case one launches multiple instances of the same
	   ! program in parallel.
	   CALL SYSTEM_CLOCK(count)
	   IF (count /= 0) THEN
		  t = TRANSFER(count, t)
	   ELSE
		  CALL DATE_AND_TIME(values=dt)
		  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
			   + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
			   + dt(3) * 24 * 60 * 60 * 60 * 1000 &
			   + dt(5) * 60 * 60 * 1000 &
			   + dt(6) * 60 * 1000 + dt(7) * 1000 &
			   + dt(8)
		  t = TRANSFER(tms, t)
	   END IF
	   s = IEOR(t(1), t(2))
	   pid = getpid() + 1099279 ! Add a prime
	   s = IEOR(s, pid)
	   IF (n >= 3) THEN
		  seed(1) = t(1) + 36269
		  seed(2) = t(2) + 72551
		  seed(3) = pid
		  IF (n > 3) THEN
			 seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
		  END IF
	   ELSE
		  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
	   END IF
	END IF
	CALL RANDOM_SEED(put=seed)
  END SUBROUTINE init_random_seed
END PROGRAM
