MODULE subs 
!
USE prec
USE global
USE inout
IMPLICIT NONE
!
CONTAINS
!
SUBROUTINE ghost
!
!       This subroutines sets up the ghost cells
!         at the ends for periodic BCs
!
  buffer = 3                 ! Will allow for 2 ghost cells on each end
  N = Nx + buffer + buffer   ! tot nodes = physical nodes + ghost nodes
  fstn = 1 + buffer          ! the first node
  lstn = Nx + buffer         ! the last node
!
  fstj = 1 + buffer          ! the first jun
  lstj = Nx + 1 + buffer     ! the last jun
!
END SUBROUTINE ghost
SUBROUTINE grid1d
!
!       Arguements:
!
!
!       Local Variables
!
  REAL(dp) :: L                                  ! Length of grid
  INTEGER :: i                                   ! Row index
  INTEGER :: j                                   ! Col index
!
!
!       Force xs < xe
!
  IF (xs > xe) THEN
    dx = xe
    xe = xs
    xs = dx
  ENDIF

  L = xe - xs
  dx = L / Nx
!
!       Create grid
!
  xj(fstj) = xs
  xc(fstn) = xs + 0.5*dx
!
  DO i = fstn+1, lstn
    xc(i) = xc(i-1) + dx
    xj(i) = xj(i-1) + dx
  ENDDO
!
    xj(lstj) = xj(lstj - 1) + dx
!
  DO i = 1, buffer
    xc(fstn - i) = xc(fstn) - real(i)*dx
    xc(lstn + i) = xc(lstn) + real(i)*dx
!
    xj(fstj - i) = xj(fstj) - real(i)*dx
    xj(lstj + i) = xj(lstj) + real(i)*dx
  END DO
!
!
END SUBROUTINE grid1d
SUBROUTINE init
!
!       This subroutine will put the initial condition on the
!        physical domain
!
!       Arguements:
!
!!!  REAL(dp), INTENT(IN), DIMENSION(N) :: x           ! Grid
!!!  REAL(dp), INTENT(OUT), DIMENSION(N) :: u          ! Function
!!!  LOGICAL, INTENT(OUT) :: ok                     ! True on success
!
!       Local Variables
!
  INTEGER :: i, j
  REAL(dp), DIMENSION(N) :: dwork                   ! soliton solutions
  REAL(dp) :: g, xwv
!
!
!
    ao = 0.0
    uo = 0.0
!
  SELECT CASE (ic)
!
    CASE (1)
!
!       MMS 1
!
      DO i = fstn, lstn
        ao(i) = 0.5 + 0.25*SIN(xc(i) - tstart)
      END DO
        uo = 2.0
!
    CASE (2)
!
!       MMS 2
!
        ao = 0.5
      DO j = fstj, lstj
        uo(j) = 2.0 - SIN(xj(j) - tstart)
      END DO
!
    CASE (3)
!
!       MMS 3
!
      DO i = fstn, lstn
        ao(i) = 10.0 + SIN(xc(i) - tstart)
      END DO
!
      DO j = fstj, lstj
        uo(j) = 1.0 - 0.5*SIN(xj(j) - tstart)
      END DO
!
    CASE (33)
!
!       MMS 33 (aka 3b)
!
      DO i = fstn, lstn
!swap        ao(i) = 3.0 - SIN(xc(i) - tstart)
        ao(i) = 0.5 + 0.25*SIN(xc(i) - tstart)
      END DO
!
      DO j = fstj, lstj
!swap        uo(j) = 0.5 + 0.25*SIN(xj(j) - tstart)
        uo(j) = 3.0 - SIN(xj(j) - tstart)
      END DO
!
    CASE (4)
!
!       MMS 4
!
      g = 10.0*(1. - EXP(-time))
!
      DO i = 1, N
        ao(i) = 8.0/SQRT(100. + 2.*g*xc(i))
      END DO
!
      DO j = 1, N+1
        uo(j) = SQRT(100. + 2.*g*xj(j))
      END DO
!
    CASE (101)
!
!       Water faucet
!
      ao = 0.8
      uo = 10.0
!
    CASE (102)
!
!       Dam break
!
      ao = 0.1
      uo = 0.0
      DO i = 1, buffer + INT(Nx/2)
        ao(i) = 10.
      END DO
!
    CASE (201)
!
!       Kreiss Problem
!
      DO i = 1, N
        ao(i) = EXP(-2. *xc(i)**2)
      END DO
!
      DO j = 1, N+1
        uo(j) = EXP(-4. *xj(j)**2)
      END DO
!
    CASE (202)
!
!       Holmas Problem
!
      uo = 1.0
      DO i = 1, N
        xwv = 0.0
        IF (ABS(xc(i) - 0.35) .lt. 0.05) xwv = 1.0
        ao(i) = 0.5 &
              + xwv*0.05*SIN(2.*pi_d*(xc(i) - 0.3)/0.1) &
                 + 0.005*SIN(2.*pi_d*xc(i)/0.01)
      END DO
    CASE (203)
!
!
!
      uo = 0.1
      ao = 0.5
      DO i = 1, N
        IF (ABS(xc(i) - 0.25) .lt. 0.05) THEN
          ao(i) = ao(i) + 0.05*SIN(2.*pi_d*(xc(i) - 0.2)/0.1)
!!!        ao(i) = 0.5 + 0.01*EXP(-2.*xc(i)**2)
        END IF
      END DO
    CASE (204)
!
!       Kinematic Wave
!
      DO i = 1, N
        ao(i) = EXP(-2.*xc(i)**2)
      END DO
!
      DO j = 1, N+1
        uo(j) = 2.0 + EXP(-4.*xj(j)**2)
      END DO
!
    CASE DEFAULT
      ao = 0.0
      uo = 0.0
  END SELECT
!
!
END SUBROUTINE init
SUBROUTINE update
!
!       Update the variables
!
  INTEGER :: c, j
!
  DO c = 1, N
    ao(c)     = an(c)
  END DO
!
! 
  DO j = 1, N+1
    uo(j)     = un(j)
  END DO
!
END SUBROUTINE update
SUBROUTINE bndry(a,u)
!
!       Copies ends of physical domain into the
!        ghost cells on the ends
!
!!!   REAL(dp), INTENT(IN) :: t                     ! time
  REAL(dp), INTENT(INOUT), DIMENSION(N) :: a  ! solution variable a
  REAL(dp), INTENT(INOUT), DIMENSION(N+1) :: u    ! solution variable u
!
!       Local variables
!
  INTEGER i, j, k
  REAL(dp) :: g
!
!       Move end of physical domain to beginning ghost cells
!
  SELECT CASE (bc)
!
!       Periodic
!
    CASE (1)
      DO k = 1, buffer
        i = k
        a(fstn - i) = a(lstn - i + 1)
        a(lstn + i) = a(fstn + i - 1)
!
        j = k
        u(fstj - j) = u(lstj - j)
        u(lstj + j) = u(fstj + j)
      END DO
!
!       Reflective
!
    CASE (2)
      DO k = 1, buffer
        i = k
        a(fstn - i) = a(fstn + i - 1)
        a(lstn + i) = a(lstn - i + 1)
!
        j = k
        u(fstj - j) = u(fstj + j)
        u(lstj + j) = u(lstj - j)
      END DO
!
!       Homogeneous Dirichlet
!
    CASE (3)
      u(fstj) = 0.0
      u(lstj) = 0.0
      DO k = 1, buffer
        i = k
        a(fstn - i) = 0.0
        a(lstn + i) = 0.0
!
        j = k
        u(fstj - j) = 0.0
        u(lstj + j) = 0.0
      END DO
!
!       Water Faucet
!
    CASE (4)
!
!       Problem dependant :: breakdown below
!
      SELECT CASE (ic)
!
!       MMS 4
!
!
        CASE (4)
!
          u(fstj) = 10.0
!
          g = 10.0*(1. - EXP(-time))
!
          DO k = 1, buffer
           i = k
           a(fstn - i) = 8.0/SQRT(100. + 2.*g*xc(fstn - i))
           a(lstn + i) = 8.0/SQRT(100. + 2.*g*xc(lstn + i))
!
           j = k
           u(fstj - j) = SQRT(100. + 2.*g*xj(fstj - j))
           u(lstj + j) = SQRT(100. + 2.*g*xj(lstj + j))
          END DO
!
        CASE (101)
!
!       Water faucet
!
          u(fstj) = 10.0
          DO k = 1, buffer
           i = k
           a(fstn - i) = 0.8
           a(lstn + i) = a(lstn)
!
           j = k
!
           u(fstj - j) = 10.0
           u(lstj + j) = u(lstj)
          END DO
!
      END SELECT    ! ic (bc = 4)
!
!
  END SELECT  ! bc
!
!
END SUBROUTINE bndry
SUBROUTINE exact
!
  INTEGER :: i, j
  REAL(dp) :: xstar, xprime, g
!
!       For special cases where an exact solution exits
!
  aex = 0.0
  uex = 0.0
!
  SELECT CASE (ic)
   CASE(1)
!
! MMS 1
!
    DO i = fstn, lstn
      aex(i) = 0.5 + 0.25*SIN(xc(i) - time)
    END DO
      uex = 2.0
!!!      CALL bndry(aex,uex)
!
   CASE (2)
!
!   MMS 2
!
      aex = 0.5
    DO j = fstj, lstj
      uex(j) = 2.0 - SIN(xj(j) - time)
    END DO
!!!      CALL bndry(aex,uex)
!
   CASE (3)
!
!   MMS 3
!
    DO i = fstn, lstn
      aex(i) = 10 + SIN(xc(i) - time)
    END DO
!
    DO j = fstj, lstj
      uex(j) = 1.0 - 0.5*SIN(xj(j) - time)
    END DO
!
!!!    CALL bndry(aex,uex)
!
!
   CASE (33)
!
!   MMS 33 (aka 3b)
!
    DO i = fstn, lstn
      aex(i) = 0.5 + 0.25*SIN(xc(i) - time)
!swap      aex(i) = 3.0 - SIN(xc(i) - time)
    END DO
!
    DO j = fstj, lstj
!swap      uex(j) = 0.5 + 0.25*SIN(xj(j) - time)
      uex(j) = 3.0 - SIN(xj(j) - time)
    END DO
!
!!!    CALL bndry(aex,uex)
!
   CASE (4)
!
!   MMS 4
!
    g = 10.0*(1. - EXP(-time))
!
    DO i = 1, N
      aex(i) = 8.0/SQRT(100. + 2.*g*xc(i))
    END DO
!
    DO j = 1, N+1
      uex(j) = SQRT(100. + 2.*g*xj(j))
    END DO
!
!
   CASE (101)
!
!  Water faucet problem
!
    xstar = xs + 10.0*time + 0.5*9.8*time**2
    DO i = fstn, lstn
      IF (xc(i) .LT. xstar) THEN
        aex(i) = 10.0*0.8/SQRT(10.0**2 + 2.*9.8*(xc(i) - xs))
      ELSE
        aex(i) = 0.8
      END IF
    END DO
    DO j = fstj, lstj
      xprime = MIN(xj(j),xstar)
      uex(j) = SQRT(10.0**2 + 2.*9.8*(xprime - xs))
    END DO
!
!
!
  END SELECT
!
!
!
END SUBROUTINE exact
SUBROUTINE source(tstar)
!
!       This subroutine updates the sources Sa, Su
!
!       Arguments:
!
  REAL(dp), INTENT(IN) :: tstar
!
!       Local
!
  INTEGER :: i, j
  REAL(dp) :: uhat, theta
  REAL(dp) :: aave
!
!
!       CASES: 0 = WF, 1 = MMS1, 2 = MMS2, 3 = MMS3
!
!
  SELECT CASE (ic)
!!!   CASE (0)
!!!!
!!!!  water faucet
!!!!
!!!    Sa = 0.0
!!!    Su = 9.8
!
   CASE (1)
!
!  MMS 1
!
    DO i = fstn, lstn
      Sa(i) = 0.25*COS(xc(i) - tstar) + 0.25*eps*SIN(xc(i) - tstar)
    END DO
!
    DO j = fstj, lstj
      Su(j) = 0.25*(sigma - hyd)*COS(xj(j) - tstar)
    END DO
!
   CASE (2)
!
!  MMS 2
!
    DO i = fstn, lstn
      Sa(i) = -0.5*COS(xc(i) - tstar)
    END DO
!
    DO j = fstj, lstj
      Su(j) = COS(xj(j) - tstar)*(SIN(xj(j) - tstar) - 1.0) &
            - mu*SIN(xj(j) - tstar)
    END DO
!
   CASE (3)
!
!  MMS 3
!
    DO i = fstn, lstn
      Sa(i) = -5.0*COS(xc(i) - tstar) &
            + eps*SIN(xc(i) - tstar) &
            - 0.5*SIN(2.0*(xc(i) - tstar))
    END DO
!
    DO j = fstj, lstj
      Su(j) = (sigma - hyd)*COS(xj(j) - tstar) &
            - 0.5*mu*SIN(xj(j) - tstar) &
            + 0.125*SIN(2.0*(xj(j) - tstar))
    END DO
!
   CASE (99)
!
!  99 = SWAP swap
!
    DO i = fstn, lstn
      Sa(i) = COS(xc(i) - tstar) &
            + 0.25*eps*SIN(xc(i) - tstar) &
            + 0.25*SIN(2.0*(xc(i) - tstar))
    END DO
!
    DO j = fstj, lstj
      Su(j) = 0.25*(8.0 - hyd + sigma)*COS(xj(j) - tstar) &
            + mu*SIN(xj(j) - tstar) &
            + 0.5*SIN(2.0*(xj(j) - tstar))    
    END DO
!
   CASE (33)
!
!  MMS 33 (aka 3b)
!
    DO i = fstn, lstn
      Sa(i) = eps*0.25*SIN(xc(i) - tstar) &
            - 0.25*SIN(2.0*(xc(i) - tstar))
    END DO
!
    DO j = fstj, lstj
      Su(j) = 0.25*(sigma - hyd - 8.0)*COS(xj(j) - tstar) &
            - 1.0*mu*SIN(xj(j) - tstar) &
            + 0.5*SIN(2.0*(xj(j) - tstar)) & 
            + fi*(3.0 - SIN(xj(j) - tstar))**2 / (0.5 &
                              - 0.25*SIN(xj(j) - tstar) )
    END DO
!
   CASE (4)
!
    DO i = fstn, lstn
      uhat = SQRT(100. + 20.*xc(i)*(1. - EXP(-tstar)))
!
      Sa(i) = -80.*EXP(-tstar)*xc(i)/(uhat**3)
    END DO
!
    DO j = fstj, lstj
      uhat = SQRT(100. + 20.*xj(j)*(1. - EXP(-tstar)))
!
      Su(j) = 10.*(1. - EXP(-tstar)) &
            + 10.*EXP(-tstar)*xj(j)/uhat
    END DO
!
   CASE (101)
!
!  water faucet
!
    Sa = 0.0
!!!    Su = 9.8*(1. - EXP(-460.*tstar))  !! Was to test similarity with MMS 4
    Su = 9.8
!
   CASE (204)
!
!  Dr B's Kinematic Wave
!
   DO i = fstn, lstn
     Sa(i) = 0.0 
   END DO
!
   DO j = fstj, lstj
     Su(j) = 0.8 
   END DO
!
   CASE DEFAULT
!
!  null
!
    Sa = 0.0
    Su = 0.0
  END SELECT
!
!       for periodic bc's wrap the sources too
!
!!!  if (bc .EQ. 1)  CALL bndry(Sa,Su)   ! in main now
!
!
END SUBROUTINE source
SUBROUTINE fou
!
!
!   Donors the void fraction to the junctions with FOU
!
  INTEGER :: i, j, Lj, Rj
!
  DO j = fstj, lstj
!
!       indices
!
    i = j
    Lj = i - 1
    Rj = i
!
    IF (uo(j) .LT. 0.0) THEN
      ajd(j) = an(Rj)
    ELSE
      ajd(j) = an(Lj)
    END IF
!
  END DO
!
END SUBROUTINE fou
SUBROUTINE mma
!
!
!   Donors the void fraction to the junctions with the GPL flux limiter
!
  INTEGER :: c, j, LLc, Lc, Rc, RRc
  REAL(dp) :: dau, daf, r, f, small
!
!
! index:
!
!       |       |       |       |       |
!       |  c-2  |  c-1  |   c   |  c+1  |  c+2
!       |       |       |       |       |
!
!      j-2     j-1      j      j+1     j+2
!
  small = 1.0d-6
!
  DO j = fstj, lstj
!
    c = j
!
    LLc = c - 2
    Lc  = c - 1
    Rc  = c
    RRc = c + 1
!
    IF (un(j) .GT. 0.0) THEN
!
      daf = an(Rc) - an(Lc)
      dau = an(Lc) - an(LLc)
!
      IF (ABS(dau) .LT. small) THEN
!
!               U > 0
!
        IF (daf .GT. SMALL) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = daf / dau
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ajd(j) = an(Lc) + 0.5*f*(an(Lc) - an(LLc))
!!
!!
    ELSEIF (un(j) .LT. 0.0) THEN
!
!               U < 0
!
      daf = an(Rc) - an(Lc)
      dau = an(RRc) - an(Rc)
!
      IF (ABS(dau) .LT. small) THEN
        IF (daf .GT. small) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = daf / dau
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ajd(j) = an(Rc) - 0.5*f*(an(RRc) - an(Rc))
!!
!!
    ELSE
!
!       U == 0
!
      ajd(j) = 0.5*(an(Lc) + an(Rc))
!
!
    END IF
!
!
  END DO
!
END SUBROUTINE mma
SUBROUTINE mmu 
!
!
!
!
  INTEGER :: c, j, LLj, Lj, Rj, RRj
  REAL(dp) :: duu, duf, r, f, small
  REAL(dp) :: uave
!
!
!  scheme:    minmod
!    fla =
!    flm =
!
!
! index:
!
!       |       |       |       |       |
!       |  c-2  |  c-1  |   c   |  c+1  |  c+2
!       |       |       |       |       |
!
!      j-2     j-1      j      j+1     j+2
!
  small = 1.0d-6
!
  DO c = fstn-1, lstn+1
!
    j = c
!
    LLj = j - 1
    Lj  = j
    Rj  = j + 1
    RRj = j + 2
!
!!!    uave = un(Rj) * un(Lj)
!!!    IF (uave .GT. 0.0) THEN
!!!      uave = 0.5*(un(Rj) + un(Lj))
!!!    ELSE
!!!!!!      WRITE(*,*) 'Logic for u about 0 is poor'
!!!      uave = 0.0
!!!    END IF
!
    uave = 0.5*(un(Rj) + un(Lj))
!
    IF (uave .GT. 0.0) THEN
!
      duf = un(Rj) - un(Lj)
      duu = un(Lj) - un(LLj)
!
      IF (ABS(duu) .LT. small) THEN
!
!               U > 0
!
        IF (duf .GT. SMALL) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = duf / duu
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ucd(c) = un(Lj) + 0.5*f*(un(Lj) - un(LLj))
!!
!!
    ELSEIF (uave .LT. 0.0) THEN
!
!               U < 0
!
      duf = un(Rj) - un(Lj)
      duu = un(RRj) - un(Rj)
!
      IF (ABS(duu) .LT. small) THEN
        IF (duf .GT. small) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = duf / duu
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ucd(c) = un(Rj) - 0.5*f*(un(RRj) - un(Rj))
!!
!!
    ELSE
!
!       U == 0
!
      ucd(c) = 0.5*(un(Lj) + un(Rj))
!
!
    END IF
!
!
  END DO
!
END SUBROUTINE mmu
SUBROUTINE edfou2
!
!   E        D         F     O     U
!   Explicit Decoupled First Order Upwind
!
!   2 = vel first, *new vel in void
!
!
!       Arguements:
!
!
!
  INTEGER ::  i, j, Li, Ri, LLj, Lj, Rj, RRj
  REAL(dp) :: conv
!
!
!       Eq 1 "vel" update ::
!
!
  DO j = fstj + bcL, lstj - bcR
!
!       indices:
!
    i = j
    LLj = i - 2
    Lj  = i - 1
    Rj  = i
    RRj = i + 1
!
!
!!      conv = uo(j)*(ucd(Rj) - ucd(Lj))/dx
!
   conv = 0.5*(uo(j) - ABS(uo(j)))*(uo(j+1) - uo(j))/dx + &
           0.5*(uo(j) + ABS(uo(j)))*(uo(j) - uo(j-1))/dx
!!!!
!!!!    conv = 0.5*(uo(j)**2 - uo(j-1)**2)/dx
!
    un(j) = uo(j) + dt * (-conv  &
          + (hyd/dx)*(ao(Rj) - ao(Lj)) &
          + (mu/dx**2)*(uo(j+1) - 2.*uo(j) + uo(j-1)) &
          + (sigma/dx**3)*(ao(RRj) - 3.*ao(Rj) + 3.*ao(Lj) - ao(LLj)) &
          - fi*uo(j)**2/ajd(j) &
          + Su(j) )
  END DO
!
!
!  Eq 1 "void" update ::
!
!
  DO i = fstn, lstn
!
! indices
!
    j = i      ! jun index
    Li = j     ! Left jun
    Ri = j+1   ! Right jun
!
    an(i) = ao(i) + &
            dt*((-1.0/dx)*( ajd(Ri)*un(Ri) - ajd(Li)*un(Li) ) &
          + (eps/dx**2)*(ao(i+1) - 2.*ao(i) + ao(i-1)) &
          + Sa(i) &
          - kreiss*((1.0 - ao(i)/2.)*(un(Ri) - un(Li))/dx + 2.*ao(i)) &
          )
  END DO
!
END SUBROUTINE edfou2
SUBROUTINE getrhs
!
!       Get the RHS's of alpha and u equations with input values (a,u)
!
!       Arguements:
!
!!!  REAL(dp), INTENT(IN), DIMENSION(N)   :: ac    ! *current* a value
!!!  REAL(dp), INTENT(IN), DIMENSION(N+1) :: uc    ! *current* u value
!
!       Local:
!
  INTEGER ::  c, j, Lj, Rj, LLc, Lc, Rc, RRc
  REAL(dp) :: conv
!
!
!  Eq 1 "void" update ::
!
!
  DO c = fstn, lstn
!
! indices
!
    j = c      ! jun index
    Lj = j     ! Left jun
    Rj = j+1   ! Right jun
!
    arhs(c) = ((-1.0/dx)*( ajd(Rj)*un(Rj) - ajd(Lj)*un(Lj) ) &
            + (eps/dx**2)*(an(c+1) - 2.*an(c) + an(c-1)) &
            + Sa(c) &
            - kreiss*((1.0 - an(c)/2.)*(un(Rj) - un(Lj))/dx + 2.*an(c)) &
            )
  END DO
!
!
!       Eq 2 "vel" update ::
!
!
  DO j = fstj + bcL, lstj - bcR
!
!       indices:
!
    c = j
    LLc = c - 2
    Lc  = c - 1
    Rc  = c
    RRc = c + 1
!
!
      conv = un(j)*(ucd(Rc) - ucd(Lc))/dx
!
!!   conv = 0.5*(un(j) - ABS(un(j)))*(un(j+1) - un(j))/dx + &
!!           0.5*(un(j) + ABS(un(j)))*(un(j) - un(j-1))/dx
!
   urhs(j) = (-conv  &
           + (hyd/dx)*(an(Rc) - an(Lc)) &
           + (mu/dx**2)*(un(j+1) - 2.*un(j) + un(j-1)) &
           + (sigma/dx**3)*(an(RRc) - 3.*an(Rc) + 3.*an(Lc) - an(LLc)) &
           - fi*un(j)**2/ajd(j) &
           + Su(j) )
  END DO
!
!
!
!
!
END SUBROUTINE getrhs
SUBROUTINE rk3s1 
!
!
!       SSPRK(3-3) Stage 1
!
!
  INTEGER :: c, j 
!
!
  DO c = fstn, lstn
    a1(c) = ao(c) + dt*arhs(c)
  END DO
!
!
  DO j = fstj, lstj
    u1(j) = uo(j) + dt*urhs(j) 
  END DO
!
!
END SUBROUTINE rk3s1
SUBROUTINE rk3s2 
!
!
!       SSPRK(3-3) Stage 2
!
!
  INTEGER :: c, j 
!
!
  DO c = fstn, lstn
    a2(c) = (3.0*ao(c) + a1(c) + dt*arhs(c))/4.0  
  END DO
!
!
  DO j = fstj, lstj
    u2(j) = (3.0*uo(j) + u1(j) + dt*urhs(j))/4.0 
  END DO
!
!
END SUBROUTINE rk3s2 
SUBROUTINE rk3s3 
!
!
!       SSPRK(3-3) Stage 3
!
!
  INTEGER :: c, j 
!
!
  DO c = fstn, lstn
    an(c) = (ao(c) + 2.0*a2(c) + 2.0*dt*arhs(c))/3.0 
  END DO
!
!
  DO j = fstj, lstj
    un(j) = (uo(j) + 2.0*u2(j) + 2.0*dt*urhs(j))/3.0 
  END DO
!
!
END SUBROUTINE rk3s3 
END MODULE subs 
