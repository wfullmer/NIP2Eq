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
