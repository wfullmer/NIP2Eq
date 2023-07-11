SUBROUTINE edfou2
!
!   E        D         F     O     U
!   Explicit Decoupled First Order Upwind
!
!   2 = vel first, *new vel in void
!
!
!       Arguements:
!
!
!
  INTEGER ::  i, j, Li, Ri, LLj, Lj, Rj, RRj
  REAL(dp) :: conv
!
!
!       Eq 1 "vel" update ::
!
!
  DO j = fstj + bcL, lstj - bcR
!
!       indices:
!
    i = j
    LLj = i - 2
    Lj  = i - 1
    Rj  = i
    RRj = i + 1
!
!
!!      conv = uo(j)*(ucd(Rj) - ucd(Lj))/dx
!
   conv = 0.5*(uo(j) - ABS(uo(j)))*(uo(j+1) - uo(j))/dx + &
           0.5*(uo(j) + ABS(uo(j)))*(uo(j) - uo(j-1))/dx
!!!!
!!!!    conv = 0.5*(uo(j)**2 - uo(j-1)**2)/dx
!
    un(j) = uo(j) + dt * (-conv  &
          + (hyd/dx)*(ao(Rj) - ao(Lj)) &
          + (mu/dx**2)*(uo(j+1) - 2.*uo(j) + uo(j-1)) &
          + (sigma/dx**3)*(ao(RRj) - 3.*ao(Rj) + 3.*ao(Lj) - ao(LLj)) &
          - fi*uo(j)**2/ajd(j) &
          + Su(j) )
  END DO
!
!
!  Eq 1 "void" update ::
!
!
  DO i = fstn, lstn
!
! indices
!
    j = i      ! jun index
    Li = j     ! Left jun
    Ri = j+1   ! Right jun
!
    an(i) = ao(i) + &
            dt*((-1.0/dx)*( ajd(Ri)*un(Ri) - ajd(Li)*un(Li) ) &
          + (eps/dx**2)*(ao(i+1) - 2.*ao(i) + ao(i-1)) &
          + Sa(i) &
          - kreiss*((1.0 - ao(i)/2.)*(un(Ri) - un(Li))/dx + 2.*ao(i)) &
          )
  END DO
!
END SUBROUTINE edfou2
