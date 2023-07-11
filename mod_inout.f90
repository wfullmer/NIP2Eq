MODULE inout
!
!	This module contains the input/output subroutines: 
!		readin
!		amovout
!		umovout
!		plot
!
USE prec
USE global
!
IMPLICIT NONE
!
CONTAINS
!
  SUBROUTINE readin
!
!       Prompt user for input on grid info,
!       time info, coeffs, IC and solution method
!
!       Local data
!
!
  gplflag = 0  ! used in IF whether or not it's std input
  kreiss = 0.0 ! unless ic == 201
!
  WRITE(*,*) '  '
  WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(*,*) '      THIS PROGRAM SOLVES AN ILL-POSED TWO-EQUATION SYSTEM      '
  WRITE(*,*) '             a_t + u a_x + a u_x = (eps) a_xx + Sa              '
  WRITE(*,*) '     u_t + u u_x = (c) a_x + (mu) u_xx + (sigma) a_xxx + Su     '
  WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(*,*) '  '
  WRITE(*,*) '   Enter the coefficient, c :  '
  WRITE(*,*) '     c < 0 hyperbolic '
  WRITE(*,*) '     c = 0 parabolic '
  WRITE(*,*) '     c > 0 elliptic '
  WRITE(*,*) '  '
  READ(*,*) hyd
  WRITE(*,*) '  ' 
  WRITE(*,*) '   Enter the 2rd Order Coefficient, eps : '
  WRITE(*,*) '  '
  READ(*,*) eps
  WRITE(*,*) '  '
  WRITE(*,*) '   Enter the 2rd Order Coefficient, mu :  '
  WRITE(*,*) '  '
  READ(*,*) mu
  WRITE(*,*) '   Enter the 3rd Order Coefficient, sigma : '
  WRITE(*,*) '  '
  READ(*,*) sigma
  WRITE(*,*) '   Enter the friction Coefficient, fi : ' 
  WRITE(*,*) '  '
  READ(*,*) fi 
!
  WRITE(*,*) '  '
  WRITE(*,*) '            How do you want to solve the equation ???           '
  WRITE(*,*) '  '
  WRITE(*,*) '   ENTER  1    for  Explicit FOU '
  WRITE(*,*) '   ENTER  2    for  Explicit FOU w/ S(n+1) '
  WRITE(*,*) '  '
  WRITE(*,*) '   ENTER  302  for  O(2) A-B - A-M Pre-Corr w/ GPL '
  WRITE(*,*) '   ENTER  303  for  O(3) A-B - A-M Pre-Corr w/ GPL '
  WRITE(*,*) '  '
  WRITE(*,*) '   ENTER  403  for  O(3) SSP RK w/ GPL '
  WRITE(*,*) '  '
  READ(*,*) method
  IF (method .GE. 300 .AND. method .LE. 499) THEN
    WRITE(*,*) '  '
    WRITE(*,*) '        WARNING :: ONLY IMPLEMENTED FOR PERIODIC BCs :: bc=1  '
    WRITE(*,*) '        WARNING :: LOGIC FOR u ABOUT 0 IS QUESTIONABLE        '
    WRITE(*,*) '  '
    WRITE(*,*) '     ENTER 1 for minmod   '
    WRITE(*,*) '     ENTER 2 for muscl    '
    WRITE(*,*) '     ENTER 3 for smart    '
    WRITE(*,*) '  '
    READ(*,*) gplflag
    IF (gplflag .EQ. 1) THEN
      fla = -1.0
      flk = -1.0
      flm =  1.0
    ELSEIF (gplflag .EQ. 2) THEN
      fla =  0.0
      flk =  0.0
      flm =  2.0
    ELSEIF (gplflag .EQ. 3) THEN
      fla =  0.0
      flk =  0.5
      flm =  4.0
    ELSE
      WRITE(*,*) '  '
      WRITE(*,*) '   INVALID ENTRY.... SWITCHING TO EDFOU2 '
      WRITE(*,*) '  '
      gplflag = 999
    END IF
  END IF
!
  IF (gplflag .EQ. 999) method = 2
!
!
  WRITE(*,*) '              What Problem do you want to solve     ???         '
  WRITE(*,*) '  '
  WRITE(*,*) '                          Verification                          '
  WRITE(*,*) '  '
  WRITE(*,*) '   ENTER   1     for     MMS 1 : a = SIN, u = 2  ' 
  WRITE(*,*) '   ENTER   2     for     MMS 2 : a = 0.5, u = SIN '
  WRITE(*,*) '   ENTER   3     for     MMS 3 : a = SIN, u = SIN '
  WRITE(*,*) '   ENTER   4     for     MMS 4 : a,u = f( EXP(t), x^ -/+ .5 ) ' 
  WRITE(*,*) '   ENTER  33     for     MMS 3b: variant of 3 for kinematic wave '
  WRITE(*,*) '  '
  WRITE(*,*) '                           Validation                           '
  WRITE(*,*) '  '
  WRITE(*,*) '   ENTER   101     for     Water Faucet ' 
  WRITE(*,*) '   ENTER   102     for     Dam Break  ' 
  WRITE(*,*) '  '
  WRITE(*,*) '                             Others                             '
  WRITE(*,*) '  '
  WRITE(*,*) '   ENTER   201     for     Kreiss-type Problem ' 
  WRITE(*,*) '   ENTER   202     for     Holmas-type Problem  ' 
  WRITE(*,*) '   ENTER   203     for '
  WRITE(*,*) '   ENTER   204     for     Dr Bs Kinematic Wave '
  WRITE(*,*) '  '
  WRITE(*,*) '                       Two-Eq TFM Problems                      '
  WRITE(*,*) '  '
  WRITE(*,*) '   ENTER   901     for     ICMF 13 stable case '
  WRITE(*,*) '   ENTER   902     for     ICMF 13 kinematically unstable case '
  WRITE(*,*) '   ENTER   903     for     ICMF 13 dnyamically unstable case '
  READ(*,*) ic
  IF (ic .EQ. 201) kreiss = 1.0
!
  WRITE(*,*) '         What boundary condition do you want to use ???         '
  WRITE(*,*) '   ENTER   1     for     Periodic '
  WRITE(*,*) '   ENTER   2     for     Mirror '
  WRITE(*,*) '   ENTER   3     for     Homogeneous Dirichlet '
  WRITE(*,*) '   ENTER   4     for     Special (IC based) '
  WRITE(*,*) '  '
  READ(*,*) bc
!
  WRITE(*,*) '  Enter the start point, end point and Num of cells for a 1D grid '
  WRITE(*,*) '  '
  WRITE(*,*) '!!!   To enter in units of pi, enter a NEGATIVE Num of cells   !!!'
  WRITE(*,*) '  '
  READ(*,*) xs
  READ(*,*) xe
  READ(*,*) Nx
  IF(Nx .lt. 0) THEN
    xs = xs * pi_d
    xe = xe * pi_d
    Nx = -one_d * Nx
  ENDIF
!
  WRITE(*,*) '                       Enter the time step                        '
  WRITE(*,*) '         !!! To enter in units of dx use NEGATIVE dt !!!          '
  WRITE(*,*) '  '
  READ(*,*) dt
!
  WRITE(*,*) '         Enter the start and end time of the simulation           '
  WRITE(*,*) '   !!! To enter in units of pi, enter a NEGATIVE END TIME !!!     '  
  WRITE(*,*) '  '
  READ(*,*) tstart
  READ(*,*) tend
  IF (tend .LT. 0) THEN
    tend = -1.0*tend*pi_d
  END IF
  WRITE(*,*) '           Enter the print frequency in number of dts             '
  WRITE(*,*) '  '
  READ(*,*) freq
  WRITE(*,*) '                        Enter the plot time                       '
  WRITE(*,*) '  '
  READ(*,*) plot_time
!
 END SUBROUTINE readin
!
!
!
 SUBROUTINE amovout(x,u,io)
!
!       This subroutine will produce a formatted output file
!         compatable with TECPLOT
!
!       Arguements:
!
  INTEGER, INTENT(IN) :: io 
  REAL(dp), INTENT(IN), DIMENSION(N) :: x, u          ! Function
!
!       Local vars:
!
  INTEGER :: i, j
!
  OPEN(UNIT=1, FILE='output1.dat', STATUS='UNKNOWN')
  OPEN(UNIT=5, FILE='output5.dat', STATUS='UNKNOWN')
  100 FORMAT(F12.5,F12.5)
!
  IF (io .EQ. 1) WRITE(io,*) 'VARIABLES = "x", "a" '
  IF (io .EQ. 5) WRITE(io,*) 'VARIABLES = "x", "aex" '
  WRITE(io,*) 'ZONE I=',Nx,', ZONETYPE=ORDERED,'
  WRITE(io,*) 'DATAPACKING=POINT, SOLUTIONTIME=',time
!
  DO i = fstn,lstn
      WRITE(io,100) x(i), u(i)
  END DO
!
!  CLOSE(io)  <<<------- REMEMBER TO CLOSE UNIT=1 POINTER IN THE MAIN PROGRAM FILE !!!
!
 END SUBROUTINE amovout
!
!
!
 SUBROUTINE umovout(x,u,io)
!
!       This subroutine will produce a formatted output file
!         compatable with TECPLOT
!
!       Arguements:
!
  INTEGER, INTENT(IN) :: io
  REAL(dp), INTENT(IN), DIMENSION(N+1) :: x, u         ! Function
!
!       Local vars:
!
  INTEGER :: i, j
!
  OPEN(UNIT=2, FILE='output2.dat', STATUS='UNKNOWN')
  OPEN(UNIT=6, FILE='output6.dat', STATUS='UNKNOWN')
  200 FORMAT(F12.5,F12.5)
!
  WRITE(io,*) 'VARIABLES = "x", "u"'
  WRITE(io,*) 'ZONE I=',Nx+1,', ZONETYPE=ORDERED,'
  WRITE(io,*) 'DATAPACKING=POINT, SOLUTIONTIME=',time
!
  DO i = fstj, lstj
      WRITE(io,200) x(i), u(i)
  END DO
!
!  CLOSE(io)  <<<------- REMEMBER TO CLOSE UNIT=2 POINTER IN THE MAIN PROGRAM FILE !!!
!
 END SUBROUTINE umovout
!
!
!
SUBROUTINE plot(x,u1,u2,io)
!
!       This subroutine will produce a formatted output file
!         compatable with ACGRACE
!
!       Arguements:
!
  INTEGER, INTENT(IN) :: io
  REAL(dp), INTENT(IN), DIMENSION(N) :: x, u1, u2          ! Function
!
!       Local vars:
!
  INTEGER :: i, j
!
  OPEN(UNIT=3, FILE='output3.dat', STATUS='UNKNOWN')
  OPEN(UNIT=4, FILE='output4.dat', STATUS='UNKNOWN')
  300 FORMAT(F12.5,F12.5,F12.5)
!
  IF (io .EQ. 3) THEN
    WRITE(io,*) '# time = ', time, '  ::  xc    void'
    DO i = 1, N
      WRITE(io,300) x(i), u1(i), u2(i)
    END DO
  ELSEIF (io .EQ. 4) THEN
    WRITE(io,*) '# time = ', time, '  ::  xj    vel'
    DO i = 1, N+1
      WRITE(io,300) x(i), u1(i), u2(i)
    END DO
  ENDIF
!
!  CLOSE(io)  <<<------- REMEMBER TO CLOSE UNIT=3,4 POINTER IN THE MAIN PROGRAM FILE !!!
!
!
 END SUBROUTINE plot
!
!
!
!
!
!
END MODULE inout
