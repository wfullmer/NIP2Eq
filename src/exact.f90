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
