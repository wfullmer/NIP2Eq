SUBROUTINE mma
!
!
!   Donors the void fraction to the junctions with the GPL flux limiter
!
  INTEGER :: c, j, LLc, Lc, Rc, RRc
  REAL(dp) :: dau, daf, r, f, small
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
  DO j = fstj, lstj
!
    c = j
!
    LLc = c - 2
    Lc  = c - 1
    Rc  = c
    RRc = c + 1
!
    IF (un(j) .GT. 0.0) THEN
!
      daf = an(Rc) - an(Lc)
      dau = an(Lc) - an(LLc)
!
      IF (ABS(dau) .LT. small) THEN
!
!               U > 0
!
        IF (daf .GT. SMALL) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = daf / dau
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ajd(j) = an(Lc) + 0.5*f*(an(Lc) - an(LLc))
!!
!!
    ELSEIF (un(j) .LT. 0.0) THEN
!
!               U < 0
!
      daf = an(Rc) - an(Lc)
      dau = an(RRc) - an(Rc)
!
      IF (ABS(dau) .LT. small) THEN
        IF (daf .GT. small) THEN
          r = flm
        ELSE
          r = 0.
        END IF
      ELSE
        r = daf / dau
      END IF
!
      f = max(0., min((2. + fla) * r, &
                       0.5*(1. + flk) * r + 0.5*(1. - flk), &
                       flm ) )
!
      ajd(j) = an(Rc) - 0.5*f*(an(RRc) - an(Rc))
!!
!!
    ELSE
!
!       U == 0
!
      ajd(j) = 0.5*(an(Lc) + an(Rc))
!
!
    END IF
!
!
  END DO
!
END SUBROUTINE mma
