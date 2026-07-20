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
