SUBROUTINE predictor
!
!
  INTEGER :: c, j 
  REAL(dp) :: ab0, ab1, ab2, ab3   ! 0 = n+1, 1 = n, ...
!
!
  ab0 = 0.0
  ab1 = 23. / 12.
  ab2 = 16. / 12.
  ab3 =  5. / 12.
!
!
!       The predictor for a
!
!
  DO c = fstn, lstn
!
!!    an(c) = ao(c) + dt*(arhso(c))
!!    an(c) = ao(c) + dt*(1.5*arhso(c) - 0.5*arhso2(c))
    an(c) = ao(c) + dt*(ab1*arhso(c) - ab2*arhso2(c) + ab3*arhso3(c))
!
  END DO
!
!
!       The predictor for u
!
!
  DO j = fstj, lstj
!
!!    un(j) = uo(j) + dt*(urhso(j))
!!    un(j) = uo(j) + dt*(1.5*urhso(j) - 0.5*urhso2(j))
    un(j) = uo(j) + dt*(ab1*urhso(j) - ab2*urhso2(j) + ab3*urhso3(j))
!
  END DO
!
!
END SUBROUTINE predictor
