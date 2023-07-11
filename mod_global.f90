MODULE global
!
!	This module contains all of the global data
!
USE prec
IMPLICIT NONE
!
!       Define Global data
!
INTEGER :: Nx         ! number of nodes in the grid in the x-dir
INTEGER :: N          ! Total number of cells including ghost cells
INTEGER :: buffer     ! # of ghost cells on each end
INTEGER :: fstn       ! first node
INTEGER :: fstj       ! first jun
INTEGER :: lstn       ! last node
INTEGER :: lstj       ! last jun
INTEGER :: method     ! selector for method to solve the advection eq
INTEGER :: ic         ! selector for which initial condition to use
INTEGER :: bc         ! selector for which boundary condition to use
INTEGER :: bcL        ! "left" "in"   boundary flag
INTEGER :: bcR        ! "right" "out" boundary flag
INTEGER :: gplflag    ! type of the flux limiter
!
REAL(dp) :: xs       ! starting point of domain
REAL(dp) :: xe       ! ending point of domain
REAL(dp) :: dx       ! mesh size of the x grid
REAL(dp) :: dt       ! time step size
REAL(dp) :: time     ! the current time
REAL(dp) :: tstart   ! the start time
REAL(dp) :: tend     ! end time of the simulation
REAL(dp) :: maxerr   ! maximum error at a given time level
REAL(dp) :: freq     ! print frequency
REAL(dp) :: plot_time ! plot time
REAL(dp) :: kreiss   ! extra terms in alpha Eq for Kreiss-Ystrom Equations
!
REAL(dp) :: hyd      ! coeff of void gradient
REAL(dp) :: eps      ! coeff of diffusion for alpha 
REAL(dp) :: mu       ! coeff of diffusion for u, kinematic viscosity
REAL(dp) :: sigma    ! coeff of dispersion for u, surface tension 
REAL(dp) :: fi       ! coeff of friction 
!
REAL(dp) :: fla, flm, flk            !GPL flags (general-piecewise flux limiter)
!
REAL(dp), ALLOCATABLE :: xc(:)       ! Centroid values for grid
REAL(dp), ALLOCATABLE :: xj(:)       ! junction values 
REAL(dp), ALLOCATABLE :: ao(:)       ! old-time values  : alpha
REAL(dp), ALLOCATABLE :: uo(:)       !                    u
REAL(dp), ALLOCATABLE :: a1(:)       ! RK stage 1 values: 
REAL(dp), ALLOCATABLE :: u1(:)       ! 
REAL(dp), ALLOCATABLE :: a2(:)       ! RK stage 2 values: 
REAL(dp), ALLOCATABLE :: u2(:)       !
REAL(dp), ALLOCATABLE :: an(:)       ! new-time values  :
REAL(dp), ALLOCATABLE :: un(:)       ! 
REAL(dp), ALLOCATABLE :: ajd(:)      ! donored void fraction   c -> j
REAL(dp), ALLOCATABLE :: ucd(:)      ! donored velocity        j -> c
!
REAL(dp), ALLOCATABLE :: arhs(:)     ! *current* RHS of a 
REAL(dp), ALLOCATABLE :: urhs(:)     ! *current* RHS of u 
!
REAL(dp), ALLOCATABLE :: Sa(:)       ! void Source
REAL(dp), ALLOCATABLE :: Su(:)       ! Vel  Source 
!
REAL(dp), ALLOCATABLE :: aex(:)      ! exact soln if needed
REAL(dp), ALLOCATABLE :: uex(:)      ! exact soln if needed
!
LOGICAL :: ok
!
!
!
END MODULE global
