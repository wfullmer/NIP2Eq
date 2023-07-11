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
      ajd(j) = an(Rj)
    ELSE
      ajd(j) = an(Lj)
    END IF
!
  END DO
!
END SUBROUTINE fou
