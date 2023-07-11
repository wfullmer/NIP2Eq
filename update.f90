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
