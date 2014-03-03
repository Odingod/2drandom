PROGRAM bd
  ! These are comments
  !
  ! Random Deposition with surface relaxation
  !
  ! Compile as: gfortran -o foo foo.f90
  !
  IMPLICIT NONE ! This forces us to define all variables
  INTEGER, ALLOCATABLE ::  h(:,:) ! height
  REAL(KIND(1.d0)), ALLOCATABLE ::  w(:) ! roughness(^2)
  INTEGER :: i1, i2, j, N, Nh, Nr, ii, jj, k
  CHARACTER(len=40) :: w_file, h_file
  LOGICAL :: do_sr=.TRUE.
  REAL :: r
  CALL init_random_seed()

  WRITE(*,*) 'Do SR? (false makes RD), say t or f '
  READ(*,*) do_sr
  WRITE(*,*) 'Number of sites (horizontal size of cluster) [try 100]'
  READ(*,*) N
  ALLOCATE(h(0:N+1,0:N+1))
  WRITE(*,*) 'How many particles (per column)? [try 1000]'
  READ(*,*) Nh
  WRITE(*,*) 'How many runs? [try 100]'
  READ(*,*) Nr
  ALLOCATE(w(Nh))

  IF (do_sr) THEN 
     WRITE(w_file,'(a,i0,a)') 'sr_w_average', N, '.txt'
     WRITE(h_file,'(a,i0,a)') 'sr_h_single', N, '.txt'
  ELSE
     WRITE(w_file,'(a,i0,a)') 'w_average', N, '.txt'  
     WRITE(h_file,'(a,i0,a)') 'h_single', N, '.txt'   
  END IF
  OPEN(123,file=w_file)
  OPEN(124,file=h_file)

  w=0.d0
  DO jj=1, Nr
     WRITE(*,*) 'Run ', jj, '/', Nr
     h=0  ! This contains heights.
     DO ii=1, Nh
        DO j=1, N
           DO k=1 , N
              i1=random_site()
              i2=random_site()
              IF (do_sr) THEN
                 IF (h(i1-1,i2)<h(i1,i2) .AND. h(i1+1,i2)<h(i1,i2) .AND. h(i1,i2-1)<h(i1,i2) .AND. h(i1,i2+1)<h(i1,i2)) THEN
                    CALL RANDOM_NUMBER(r)
                    IF (r>.75) THEN
                       i1=i1+1
                    ELSE IF (r>0.5) THEN
                       i1=i1-1
                    ELSE IF (r>0.25) THEN
                       i2=i2+1
                    ELSE 
                       i2=i2-1
                    END IF
                 ELSE IF (h(i1-1,i2)<h(i1,i2) .AND. h(i1+1,i2)<h(i1,i2) .AND. h(i1,i2-1)<h(i1,i2)) THEN !EI Y
                    CALL RANDOM_NUMBER(r)
                    IF (r>.66) THEN
                       i1=i1+1
                    ELSE IF (r>0.33) THEN
                       i1=i1-1
                    ELSE 
                       i2=i2-1
                    END IF
                 ELSE IF (h(i1-1,i2)<h(i1,i2) .AND. h(i1+1,i2)<h(i1,i2) .AND. h(i1,i2+1)<h(i1,i2)) THEN !EI A
                    CALL RANDOM_NUMBER(r)
                    IF (r>.66) THEN
                       i1=i1+1
                    ELSE IF (r>0.33) THEN
                       i1=i1-1
                    ELSE 
                       i2=i2+1
                    END IF 
                 ELSE IF (h(i1-1,i2)<h(i1,i2) .AND. h(i1,i2-1)<h(i1,i2) .AND. h(i1,i2+1)<h(i1,i2)) THEN !EI O
                    CALL RANDOM_NUMBER(r)
                    IF (r>.66) THEN
                       i2=i2+1
                    ELSE IF (r>0.33) THEN
                       i1=i1-1
                    ELSE 
                       i2=i2-1
                    END IF 
                 ELSE IF (h(i1+1,i2)<h(i1,i2) .AND. h(i1,i2-1)<h(i1,i2) .AND. h(i1,i2+1)<h(i1,i2)) THEN !EI V
                    CALL RANDOM_NUMBER(r)
                    IF (r>.66) THEN
                       i1=i1+1
                    ELSE IF (r>0.33) THEN
                       i2=i2+1
                    ELSE 
                       i2=i2-1
                    END IF
                 ELSE IF (h(i1+1,i2)<h(i1,i2) .AND. h(i1,i2-1)<h(i1,i2)) THEN !OA
                    CALL RANDOM_NUMBER(r)
                    IF (r>.5) THEN
                       i1=i1+1
                    ELSE 
                       i2=i2-1
                    END IF
                 ELSE IF (h(i1+1,i2)<h(i1,i2) .AND. h(i1,i2+1)<h(i1,i2)) THEN !OY
                    CALL RANDOM_NUMBER(r)
                    IF (r>.5) THEN
                       i1=i1+1
                    ELSE 
                       i2=i2+1
                    END IF
                 ELSE IF (h(i1+1,i2)<h(i1,i2) .AND. h(i1-1,i2)<h(i1,i2)) THEN !OV
                    CALL RANDOM_NUMBER(r)
                    IF (r>.5) THEN
                       i1=i1+1
                    ELSE 
                       i1=i1-1
                    END IF
                 ELSE IF (h(i1-1,i2)<h(i1,i2) .AND. h(i1,i2+1)<h(i1,i2)) THEN !YV
                    CALL RANDOM_NUMBER(r)
                    IF (r>.5) THEN
                       i1=i1-1
                    ELSE 
                       i2=i2+1
                    END IF
                 ELSE IF (h(i1,i2+1)<h(i1,i2) .AND. h(i1,i2-1)<h(i1,i2)) THEN !YA
                    CALL RANDOM_NUMBER(r)
                    IF (r>.5) THEN
                       i2=i2+1
                    ELSE 
                       i2=i2-1
                    END IF
                 ELSE IF (h(i1-1,i2)<h(i1,i2) .AND. h(i1,i2-1)<h(i1,i2)) THEN !VA
                    CALL RANDOM_NUMBER(r)
                    IF (r>.5) THEN
                       i1=i1-1
                    ELSE 
                       i2=i2-1
                    END IF
                 ELSE IF (h(i1-1,i2)<h(i1,i2)) THEN !V
                    i1=i1-1
                 ELSE IF (h(i1+1,i2)<h(i1,i2)) THEN !O
                    i1=i1+1
                 ELSE IF (h(i1,i2-1)<h(i1,i2)) THEN !A 
                    i2=i2-1
                 ELSE IF (h(i1,i2+1)<h(i1,i2)) THEN !Y
                    i2=i2+1
                 END IF
                    

                    
                    
               
                 
                 IF (i1==0) i1=N
                 IF (i1==N+1) i1=1
                 IF (i2==0) i2=N
                 IF (i2==N+1) i2=1

                 h(i1,i2)=h(i1,i2)+1 ! increase at selected column
                 IF (i1==1) h(N+1,i2)=h(1,i2) ! These take care of 
                 IF (i1==N) h(0,i2)=h(N,i2)   ! the periodic boundary conditions.
                 IF (i2==1) h(i1,N+1)=h(i1,1) ! These take care of 
                 IF (i2==N) h(i1,0)=h(i1,N)   ! the periodic boundary conditions.
              ELSE
                 Write(*,*) "Doing RD"
                 h(i1,i2)=h(i1,i2)+1 ! This is the rule of RD.
              END IF
           END DO
        END DO
        w(ii)=w(ii)+rough2(h(1:N,1:N))
     END DO
  END DO
  w=w/REAL(Nr,KIND(1.d0))
  DO ii=1, Nh
     WRITE(123,*) ii, SQRT(w(ii))
  END DO
  WRITE(124,*) h
  CLOSE(123)
  CLOSE(124)
  
  WRITE(*,*) 'Data in ', w_file
CONTAINS
  FUNCTION random_site()
    INTEGER :: random_site
    REAL(KIND(1.d0)) :: r
    CALL RANDOM_NUMBER(r)
    random_site=MIN(INT(r*N)+1,N)
  END FUNCTION random_site
  FUNCTION rough2(h)
    REAL(KIND(1.d0)) :: rough2, eh, eh2
    INTEGER :: h(:,:), N
    REAL(KIND(1.d0)), DIMENSION(SIZE(h)) :: nh
    N=SIZE(h)
    nh=RESHAPE(h, (/N/))-MINVAL(h)
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
