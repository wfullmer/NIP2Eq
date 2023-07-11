PROGRAM NIP2EQ 
!
! Purpose: 
!		NIP2EQ (quasi-Nonlinear Ill-Posed 2-equation model) will 
!	solve a (potentially) non-hyperbolic (i.e. ill-posed) IVP. The system mimics
!	the two-fluid model for void and velocity if c > 0, the 1-D shallow water 
!	equations (SWE) if c = 0, and becomes the water faucet problem if c = 0 and
! 	a gravity source is included in Eq 2. 
!
!		The system is solved in a manner similar to TFIT and 
!	other 1D TFM's: the nominal method is FOU with staggering. "Void" (a) is stored
!  	at the cell faces (xj) and the "velocity" (u) is stored at the cell centroids 
!	(xc). No special methods are needed for correction b/c "pressure" d.n.e.
! 
!
 USE prec
 USE global
 USE inout
 USE subs
!
 IMPLICIT NONE       ! Force explicit declaration of variables
!
! 	Data Dictionary:
!
 INTEGER :: ii, jj
 INTEGER :: dbug
 REAL(dp) :: print_time, crash, term
 REAL(dp) :: L2a, L2u
 REAL(dp) :: aamp, uamp 
 REAL(dp) :: fft_time
!
!
 dbug = 7
 OPEN(UNIT=dbug, FILE='info.dat', STATUS='UNKNOWN')
!
!
!	Get user supplied input data
!
!
 CALL readin
!
!
!       Create ghost cells for periodic BCs
!
!
 CALL ghost
!
!
!	Allocate memory
!
!
!
 ALLOCATE(xc(N),an(N),ao(N),Sa(N),aex(N))
 ALLOCATE(arhs(N),a1(N),a2(N))
!
 ALLOCATE(xj(N+1),un(N+1),uo(N+1),Su(N+1),uex(N+1))
 ALLOCATE(urhs(N+1),u1(N+1),u2(N+1))
!
 ALLOCATE(ajd(N+1),ucd(N))
!
!	Generate the grid
!
!
 CALL grid1d 
!
!       Boundary flags
!
 IF (bc .EQ. 4) THEN
  IF (ic .EQ. 4 .OR. ic .EQ. 101) THEN 
    bcL = 1                           !!!!! THIS IMPLEMENTATION WILL ONLY WORK FOR mu = eps = sigma = 0
    bcR = 0
  ELSE
    bcL = 0
    bcR = 0
  END IF
 END IF
!	Set up the time step and print frequency
! 	  if they are dx dependant
!
!
 IF (dt .LT. 0.) THEN
   dt = -dt * dx
 END IF
!
 IF (freq .LE. 0) THEN
   freq = (((tstart - tend)/dt)/freq)
 END IF
 plot_time = plot_time + dt/100.
 fft_time  = 3.1416 
!
!	Initialize
!
!
!
 L2a = 0.0
 L2u = 0.0
 ok = .TRUE.
 time = tstart
 print_time = tstart + (freq-1)*dt
!
 CALL init
 CALL bndry(ao,uo)
 an = ao
 un = uo
 a1 = ao
 u1 = uo
 a2 = ao
 u2 = uo
 ajd  = 0.0
 ucd  = 0.0
 aex  = 0.0
 uex  = 0.0
 urhs = 0.0
 arhs = 0.0
 CALL exact
 IF (ic .LE. 3) CALL bndry(aex,uex)
 IF (ic .EQ. 33) CALL bndry(aex,uex)
!
 term = 1000.*ABS(MAXVAL(an))
!
 CALL amovout(xc,an,1)
 CALL umovout(xj,un,2)
 CALL amovout(xc,aex,5)
 CALL umovout(xj,uex,6)
 WRITE(*,*) 'Movie data generated for time =', time
!
!
!       Time Iterations of NIP2Eq
!
!
 DO WHILE (time .LT. tend)
!
!
!       Update previous variables : n+1 -> n, n -> n-1, ...
!       ***CURRENT*** TIME IS n+1 level
!
!
   time = time + dt
   CALL update 
!!!   CALL bndry(ao,uo)                   ! shouldn't be needed, jic
!!!   CALL source(time) 
!!!   IF(bc .LE. 2) CALL bndry(Sa,Su)     ! shouldn't be needed, jic
!
!
!	Get new time variables:
!
!
!
  SELECT CASE (method)
   CASE (1)
!
!
!       edfou2 with old-time sources
!
!
     CALL source(time - dt)
     IF (bc .LE. 2) CALL bndry(Sa, Su)
     CALL fou
     CALL edfou2
!
   CASE (2)
!
!
!       edfou2 with new time sources
!
!
     CALL source(time)
     IF (bc .LE. 2) CALL bndry(Sa,Su)
     CALL fou
     CALL edfou2
!
   CASE (303)
!
!
!       3rd O Adams-Bashforth - Adams-Moulton Predictor-Corrector 
!       with GPL flux limiter
!
!
!!     CALL source(time - dt)   ! get sources values with old-time data
!!     CALL bndry(Sa,Su)
!!     CALL mma                 ! get donored values with old-time data
!!     CALL mmu
!!     CALL getrhsao            ! get rhs's          with old-time data
!!     CALL getrhsuo
!!     CALL predictor           ! get new values     with old-time data
!!     CALL bndry(an,un)
!!!
!!     CALL source(time)        ! get sources values with *new*-time data
!!     CALL bndry(Sa,Su)
!!     CALL mma                 ! get donored values with *new*-time data
!!     CALL mmu
!!     CALL getrhsac            ! get rhs's          with *new*-time data
!!     CALL getrhsuc            !  stored in _rhso3
!!     CALL corrector           ! get new values     with *new*-time data
!
   CASE (403)
!
!
!       3rd O SSP Runge-Kutta (Shu, Shu + Osher, ...) with GPL flux-lim
!       !!! NOTE !!!  mma, mmu, getrhs use (an,un) so store *current* variable
!       in (an,un) before updateing 
!
!
!       stage 1: _1 & f(_o,to)       
!
!    _n = _o  --> already stored here from update
     CALL source(time - dt)
     IF (bc .LE. 2) CALL bndry(Sa,Su)
     CALL mma
     CALL mmu
     CALL getrhs
     CALL rk3s1
     CALL bndry(a1,u1)
!
!       stage 2 : _2 & f(_1,tn)
!
     an = a1
     un = u1
     CALL source(time)
     IF (bc .LE. 2) CALL bndry(Sa,Su)
     CALL mma
     CALL mmu
     CALL getrhs
     CALL rk3s2
     CALL bndry(a2,u2)
!
!       stage 3 : _n & f(_2,tn-1/2)
!
     an = a2
     un = u2
     CALL source(time - 0.5*dt)
     IF (bc .LE. 2) CALL bndry(Sa,Su)
     CALL mma
     CALL mmu
     CALL getrhs
     CALL rk3s3
!!   CALL bndry(an,un)  ! below 
!
   CASE DEFAULT
!
!       default to case 2
!
     CALL source(time)
     IF (bc .LE. 2) CALL bndry(Sa, Su)
     CALL fou
     CALL edfou2
!
  END SELECT 
!
!
!       Place BC's on new variables
!
!
   CALL bndry(an,un)
!
!
!       DeBug
!
!
  IF (time .GT. fft_time) THEN
!!!    aamp = MAXVAL(an(fstn:lstn)) - MINVAL(an(fstn:lstn))
!!!    uamp = MAXVAL(un(fstj:lstj)) - MINVAL(un(fstj:lstj))
!!!    WRITE(dbug,*) time, aamp, uamp
    fft_time = fft_time + pi_d/2.0
    DO ii = fstn, lstn
      WRITE(dbug,*) an(ii)
    END DO
  END IF
!
!
!	ERROR Check & Quit if CRASH
!
!
  crash = ABS(MAXVAL(an))
  IF (crash .GT. term) THEN
    WRITE(*,*) '  '
    WRITE(*,*) 'Function value exceded MAXIMUM, terminating...'
    WRITE(*,*) '  '
    time = 2.*tend
  END IF
!
!
!	Print new variables to outputs
!
!       -> Moive
!
  IF (time .GT. print_time) THEN
    CALL exact
    CALL amovout(xc,an,1)
    CALL umovout(xj,un,2)
    CALL amovout(xc,aex,5)
    CALL umovout(xj,uex,6)
!
!!!!!
!!!!
!!!
!!    L2a = 0.0
!!    L2u = 0.0
!!    DO ii = fstn, lstn
!!      L2a = L2a + (an(ii) - aex(ii))**2
!!    END DO
!!    L2a = SQRT(L2a/Nx)
!!!
!!    DO jj = fstj, lstj
!!      L2u = L2u + (un(jj) - uex(jj))**2
!!    END DO
!!    L2u = SQRT(L2u/Nx)
!!    WRITE(dbug,*) time, L2a, L2u
!!!
!!!!
!!!!!
!
    print_time = print_time + freq*dt
    WRITE(*,*) 'Movie data generated for time =', time
  END IF
!
!       -> Plot
!
  IF (ABS(time - plot_time) .LE. dt/2.) THEN
!!!    IF (ic .LE. 3) THEN
      L2a = 0.0
      L2u = 0.0
      CALL exact
      WRITE(*,*) 'TIME = ', time
      DO ii = fstn, lstn
        L2a = L2a + (an(ii) - aex(ii))**2
!!!        L1a = L1a + ABS(an(ii) - aex(ii))
      END DO
!!!      L1a = L2a/Nx
      L2a = SQRT(L2a/Nx)
!
      DO jj = fstj, lstj
        L2u = L2u + (un(jj) - uex(jj))**2
      END DO
      L2u = SQRT(L2u/Nx)
      WRITE(*,*) '  '
      WRITE(*,*) 'L2a = ', L2a
      WRITE(*,*) 'L2u = ', L2u
      WRITE(*,*) '  '
!!!    END IF
    CALL plot(xc,an,aex,3)
    CALL plot(xj,un,uex,4)
    WRITE(*,*) '  '
    WRITE(*,*) 'Plot data generated for time =', time
    WRITE(*,*) '  '
  END IF
!
!
!       Repeat until tend
!
!
END DO
!
!
!
IF (time .GT. 1.5*tend) then
  CALL plot(xc,an,ao,3)
  CALL plot(xj,un,uo,4)
END IF
!
!
time = time - dt
!
!	Closing arguements
!
! 	1 = a mov, 2 = u mov, 3 = a plot, 4 = u plot, 5 = aex mov, 6 = uex mov
!
WRITE(*,*) '   '
WRITE(*,*) '   '
WRITE(*,*) 'output1.dat = alpha movie'
WRITE(*,*) 'output2.dat =   u   movie'
WRITE(*,*) 'output3.dat = alpha plot '
WRITE(*,*) 'output4.dat =   u   plot '
WRITE(*,*) 'output5.dat = a ex  movie'
WRITE(*,*) 'output6.dat = u ex  movie'
WRITE(*,*) '   '
!
!
WRITE(*,*) '   '
WRITE(*,*) ' DISCOUNT DOUBLE CHECK '
WRITE(*,*) ' hyd   = ', hyd
WRITE(*,*) ' eps   = ', eps
WRITE(*,*) ' mu    = ', mu
WRITE(*,*) ' sigma = ', sigma
WRITE(*,*) ' fi    = ', fi
!
!
CLOSE(1)
CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)
CLOSE(6)
CLOSE(dbug)
!
DEALLOCATE(xc,an,ao,Sa,aex)
DEALLOCATE(arhs,a1,a2)
DEALLOCATE(xj,un,uo,Su,uex)
DEALLOCATE(urhs,u1,u2)
DEALLOCATE(ajd,ucd)
!
!
END PROGRAM NIP2EQ 
