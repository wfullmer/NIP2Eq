MODULE donor
!
!       This module contains the global data and the subroutines :
!               fou
!               mm
!               mmcv
!
!
USE prec
USE global
USE subs
IMPLICIT NONE
!
!       Define Global data
!
!
CONTAINS
!
!
!
 SUBROUTINE fou
!
!
!   Donors the void fraction to the junctions with FOU
!
  INTEGER :: i, j, Lj, Rj
!
  DO j = fstj, lstj
!
!       indices
!
    i = j
    Lj = i - 1
    Rj = i
!
    IF (uo(j) .LT. 0.0) THEN
      ajd(j) = ao(Rj)
    ELSE
      ajd(j) = ao(Lj)
    END IF
!
  END DO
!
 END SUBROUTINE fou
!
!
!
 SUBROUTINE mma
!
!
!   Donors the void fraction to the junctions with the GPL flux limiter 
!
  INTEGER :: c, j, LLc, Lc, Rc, RRc
  REAL(dp) :: dau, daf, r, f, small
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
!
!
!
 SUBROUTINE mmu(pd,p)
!
!
   REAL(dp), INTENT(IN),  DIMENSION(N+1) :: p     ! flux at j
   REAL(dp), INTENT(OUT), DIMENSION(N) :: pd      ! donored flux at c
!
!       Using p and pd b/c you may want u* d/dx (u) or d/dx (.5 u*u)
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
  uave = 0.0
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
    uave = un(Lj) * un(Rj)
    IF (uave .LT. 0.) THEN
      uave = 0.0
    ELSE
      uave = 0.5*(un(Lj) + un(Rj))
    END IF
!
    IF (uave .GT. 0.0) THEN
!
      duf = p(Rj) - p(Lj)
      duu = p(Lj) - p(LLj)
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
      pd(c) = p(Lj) + 0.5*f*(p(Lj) - p(LLj))
!!
!!
    ELSEIF (uave .LT. 0.0) THEN
!
!               U < 0
!
      duf = p(Rj) - p(Lj)
      duu = p(RRj) - p(Rj)
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
      pd(c) = p(Rj) - 0.5*f*(p(RRj) - p(Rj))
!!
!!
    ELSE
!
!       U == 0
!
      pd(c) = 0.5*(p(Lj) + p(Rj))
!
!
    END IF
!
!
  END DO
!
 END SUBROUTINE mmu
!
!
!
!
!
!
END MODULE donor
