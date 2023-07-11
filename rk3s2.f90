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
