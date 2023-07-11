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
