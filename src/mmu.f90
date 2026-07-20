SUBROUTINE mmu 
!
!
!
!
  INTEGER :: c, j, LLj, Lj, Rj, RRj
  REAL(dp) :: duu, duf, r, f, small
  REAL(dp) :: uave
!
!
!  scheme:    minmod
!    fla =
!    flm =
!
!
! index:
!
!       |       |       |       |       |
!       |  c-2  |  c-1  |   c   |  c+1  |  c+2
!       |       |       |       |       |
!
!      j-2     j-1      j      j+1     j+2
!
  small = 1.0d-6
!
  DO c = fstn-1, lstn+1
!
    j = c
!
    LLj = j - 1
    Lj  = j
    Rj  = j + 1
    RRj = j + 2
!
!!!    uave = un(Rj) * un(Lj)
!!!    IF (uave .GT. 0.0) THEN
!!!      uave = 0.5*(un(Rj) + un(Lj))
!!!    ELSE
!!!!!!      WRITE(*,*) 'Logic for u about 0 is poor'
!!!      uave = 0.0
!!!    END IF
!
    uave = 0.5*(un(Rj) + un(Lj))
!
    IF (uave .GT. 0.0) THEN
!
      duf = un(Rj) - un(Lj)
      duu = un(Lj) - un(LLj)
!
      IF (ABS(duu) .LT. small) THEN
!
!               U > 0
!
        IF (duf .GT. SMALL) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = duf / duu
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ucd(c) = un(Lj) + 0.5*f*(un(Lj) - un(LLj))
!!
!!
    ELSEIF (uave .LT. 0.0) THEN
!
!               U < 0
!
      duf = un(Rj) - un(Lj)
      duu = un(RRj) - un(Rj)
!
      IF (ABS(duu) .LT. small) THEN
        IF (duf .GT. small) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = duf / duu
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ucd(c) = un(Rj) - 0.5*f*(un(RRj) - un(Rj))
!!
!!
    ELSE
!
!       U == 0
!
      ucd(c) = 0.5*(un(Lj) + un(Rj))
!
!
    END IF
!
!
  END DO
!
END SUBROUTINE mmu
