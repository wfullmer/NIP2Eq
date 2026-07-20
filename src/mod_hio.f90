MODULE hio
!
!       This module contains the global data and the subroutines :
!               abam3
!
!
USE prec
USE global
USE subs
IMPLICIT NONE
!
!       Define Global data
!
!
CONTAINS
!
!
!
 SUBROUTINE abam3 
!
!
!   Donors the void fraction to the junctions with FOU
!
  INTEGER :: i, j, Lj, Rj
!

  *** need to make sure source was called previously otherwise RHS will be wrong


  drop rhso3
  move rhso2 -> rhso3
  move rhso  -> rhso2
  
  drop rhso, rhsn

  eval RHS with previous an, un (now ao, uo)

  get intermediate new value with Adams-Bashforth with ao, uo
  (can store in yn b/c we're just going to use it for rhs-star)

  y-star = yo &
         + dt*(23./ 12.)*rhso &
         - dt*(16./ 12.)*rhso2 &
         + dt*(5. / 12.)*rhso3 


  eval rhs-star  (can store in rhso3 b/c you dont need it anymore)

  get 'true' new value with Adams-Moulton with a-star, u-star

  yn = yo &
     + dt*(5. / 12.)*rhsn-star &
     + dt*(8. / 12.)*rhso
     - dt*(1. / 12.)*rhso2
   



  DO j = fstj, lstj
!
!       indices
!
    i = j
    Lj = i - 1
    Rj = i
!
    IF (uo(j) .LT. 0.0) THEN
      ajd(j) = ao(Rj)
    ELSE
      ajd(j) = ao(Lj)
    END IF
!
  END DO
!
 END SUBROUTINE abam3
!
!
!
!
!
!
END SUBROUTINE hio
