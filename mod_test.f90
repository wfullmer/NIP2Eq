MODULE hio
!
!       This module contains the global data and the subroutines :
!               abam3
!
!
USE prec
USE global
USE subs
USE sources
USE donor
USE explicit
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
  INTEGER :: c, j
!
!       Get RHS's with old time data 
!
  CALL getarhs(ao,uo,arhso)     ! arhs only defined from fstn - lstn
  CALL geturhs(ao,uo,urhso)     ! urhs only defined from fstj - lstj
!
!       Time advancement
!
  DO c = fstn, lstn
    an(c) = ao(c) + dt*arhso(c)
  END DO
!
  DO j = fstj, lstj
    un(j) = uo(j) + dt*urhso(j)
  END DO
!
 END SUBROUTINE abam3
!
!
!
!
!
!
END MODULE hio
