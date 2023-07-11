SUBROUTINE corrector 
!
!
  INTEGER :: c, j 
  REAL(dp) :: ab0, ab1, ab2, ab3   ! 0 = n+1, 1 = n, ...
!
!
  ab0 = 5. / 12.
  ab1 = 8. / 12.
  ab2 = 1. / 12.
  ab3 = 0.0
!
!
!       The predictor for a
!
!
  DO c = fstn, lstn
!
!!    an(c) = ao(c) + dt*(0.5*arhso3(c) + 0.5*arhso(c))
    an(c) = ao(c) + dt*(ab0*arhso3(c) + ab1*arhso(c) - ab2*arhso2(c))
!
  END DO
!
!
!       The predictor for u
!
!
  DO j = fstj, lstj
!
!!    un(j) = uo(j) + dt*(0.5*urhso3(j) + 0.5*urhso(j))
    un(j) = uo(j) + dt*(ab0*urhso3(j) + ab1*urhso(j) - ab2*urhso2(j))
!
  END DO
!
!
END SUBROUTINE corrector
