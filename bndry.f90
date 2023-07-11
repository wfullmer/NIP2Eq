SUBROUTINE bndry(a,u)
!
!       Copies ends of physical domain into the
!        ghost cells on the ends
!
!!!   REAL(dp), INTENT(IN) :: t                     ! time
  REAL(dp), INTENT(INOUT), DIMENSION(N) :: a  ! solution variable a
  REAL(dp), INTENT(INOUT), DIMENSION(N+1) :: u    ! solution variable u
!
!       Local variables
!
  INTEGER i, j, k
  REAL(dp) :: g
!
!       Move end of physical domain to beginning ghost cells
!
  SELECT CASE (bc)
!
!       Periodic
!
    CASE (1)
      DO k = 1, buffer
        i = k
        a(fstn - i) = a(lstn - i + 1)
        a(lstn + i) = a(fstn + i - 1)
!
        j = k
        u(fstj - j) = u(lstj - j)
        u(lstj + j) = u(fstj + j)
      END DO
!
!       Reflective
!
    CASE (2)
      DO k = 1, buffer
        i = k
        a(fstn - i) = a(fstn + i - 1)
        a(lstn + i) = a(lstn - i + 1)
!
        j = k
        u(fstj - j) = u(fstj + j)
        u(lstj + j) = u(lstj - j)
      END DO
!
!       Homogeneous Dirichlet
!
    CASE (3)
      u(fstj) = 0.0
      u(lstj) = 0.0
      DO k = 1, buffer
        i = k
        a(fstn - i) = 0.0
        a(lstn + i) = 0.0
!
        j = k
        u(fstj - j) = 0.0
        u(lstj + j) = 0.0
      END DO
!
!       Water Faucet
!
    CASE (4)
!
!       Problem dependant :: breakdown below
!
      SELECT CASE (ic)
!
!       MMS 4
!
!
        CASE (4)
!
          u(fstj) = 10.0
!
          g = 10.0*(1. - EXP(-time))
!
          DO k = 1, buffer
           i = k
           a(fstn - i) = 8.0/SQRT(100. + 2.*g*xc(fstn - i))
           a(lstn + i) = 8.0/SQRT(100. + 2.*g*xc(lstn + i))
!
           j = k
           u(fstj - j) = SQRT(100. + 2.*g*xj(fstj - j))
           u(lstj + j) = SQRT(100. + 2.*g*xj(lstj + j))
          END DO
!
        CASE (101)
!
!       Water faucet
!
          u(fstj) = 10.0
          DO k = 1, buffer
           i = k
           a(fstn - i) = 0.8
           a(lstn + i) = a(lstn)
!
           j = k
!
           u(fstj - j) = 10.0
           u(lstj + j) = u(lstj)
          END DO
!
      END SELECT    ! ic (bc = 4)
!
!
  END SELECT  ! bc
!
!
END SUBROUTINE bndry
