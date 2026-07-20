SUBROUTINE getrhs
!
!       Get the RHS's of alpha and u equations with input values (a,u)
!
!       Arguements:
!
!!!  REAL(dp), INTENT(IN), DIMENSION(N)   :: ac    ! *current* a value
!!!  REAL(dp), INTENT(IN), DIMENSION(N+1) :: uc    ! *current* u value
!
!       Local:
!
  INTEGER ::  c, j, Lj, Rj, LLc, Lc, Rc, RRc
  REAL(dp) :: conv
!
!
!  Eq 1 "void" update ::
!
!
  DO c = fstn, lstn
!
! indices
!
    j = c      ! jun index
    Lj = j     ! Left jun
    Rj = j+1   ! Right jun
!
    arhs(c) = ((-1.0/dx)*( ajd(Rj)*un(Rj) - ajd(Lj)*un(Lj) ) &
            + (eps/dx**2)*(an(c+1) - 2.*an(c) + an(c-1)) &
            + Sa(c) &
            - kreiss*((1.0 - an(c)/2.)*(un(Rj) - un(Lj))/dx + 2.*an(c)) &
            )
  END DO
!
!
!       Eq 2 "vel" update ::
!
!
  DO j = fstj + bcL, lstj - bcR
!
!       indices:
!
    c = j
    LLc = c - 2
    Lc  = c - 1
    Rc  = c
    RRc = c + 1
!
!
      conv = un(j)*(ucd(Rc) - ucd(Lc))/dx
!
!!   conv = 0.5*(un(j) - ABS(un(j)))*(un(j+1) - un(j))/dx + &
!!           0.5*(un(j) + ABS(un(j)))*(un(j) - un(j-1))/dx
!
   urhs(j) = (-conv  &
           + (hyd/dx)*(an(Rc) - an(Lc)) &
           + (mu/dx**2)*(un(j+1) - 2.*un(j) + un(j-1)) &
           + (sigma/dx**3)*(an(RRc) - 3.*an(Rc) + 3.*an(Lc) - an(LLc)) &
           - fi*un(j)**2/ajd(j) &
           + Su(j) )
  END DO
!
!
!
!
!
END SUBROUTINE getrhs
