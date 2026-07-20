SUBROUTINE ghost
!
!       This subroutines sets up the ghost cells
!         at the ends for periodic BCs
!
  buffer = 3                 ! Will allow for 2 ghost cells on each end
  N = Nx + buffer + buffer   ! tot nodes = physical nodes + ghost nodes
  fstn = 1 + buffer          ! the first node
  lstn = Nx + buffer         ! the last node
!
  fstj = 1 + buffer          ! the first jun
  lstj = Nx + 1 + buffer     ! the last jun
!
END SUBROUTINE ghost
