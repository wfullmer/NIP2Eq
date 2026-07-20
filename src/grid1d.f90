SUBROUTINE grid1d
!
!       Arguements:
!
!
!       Local Variables
!
  REAL(dp) :: L                                  ! Length of grid
  INTEGER :: i                                   ! Row index
  INTEGER :: j                                   ! Col index
!
!
!       Force xs < xe
!
  IF (xs > xe) THEN
    dx = xe
    xe = xs
    xs = dx
  ENDIF

  L = xe - xs
  dx = L / Nx
!
!       Create grid
!
  xj(fstj) = xs
  xc(fstn) = xs + 0.5*dx
!
  DO i = fstn+1, lstn
    xc(i) = xc(i-1) + dx
    xj(i) = xj(i-1) + dx
  ENDDO
!
    xj(lstj) = xj(lstj - 1) + dx
!
  DO i = 1, buffer
    xc(fstn - i) = xc(fstn) - real(i)*dx
    xc(lstn + i) = xc(lstn) + real(i)*dx
!
    xj(fstj - i) = xj(fstj) - real(i)*dx
    xj(lstj + i) = xj(lstj) + real(i)*dx
  END DO
!
!
END SUBROUTINE grid1d
