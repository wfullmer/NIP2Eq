SUBROUTINE rk3s3 
!
!
!       SSPRK(3-3) Stage 3
!
!
  INTEGER :: c, j 
!
!
  DO c = fstn, lstn
    an(c) = (ao(c) + 2.0*a2(c) + 2.0*dt*arhs(c))/3.0 
  END DO
!
!
  DO j = fstj, lstj
    un(j) = (uo(j) + 2.0*u2(j) + 2.0*dt*urhs(j))/3.0 
  END DO
!
!
END SUBROUTINE rk3s3 
