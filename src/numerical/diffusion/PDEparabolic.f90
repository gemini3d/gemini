module PDEparabolic

!! a module for use for solving parabolic partial differential equations


!> banded and tridiagonal solvers, for now we just take everything to be banded...
use phys_consts, only: wp
use vendor_lapack95, only: gbsv!,gtsv

implicit none (type, external)

private
public :: TRBDF21D, backEuler1D


contains


function TRBDF21D(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,BCtype,dx1,dx1i)

!! SOLVE A 1D DIFFUSION PROBLEM.  IT IS EXPECTED THAT
!! GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!! THEY ARE PASSED INTO THIS ROUTINE.
!!
!! FORM OF EQUATION SOLVED IS:
!!  dT/dt + A T + B dT/dx + C d/dx(D dT/dx) = E
!!
!!  NOTE: UPON FURTHER REVIEW I THINK THE FORM SOLVED IS ACTUALLY:
!!  dT/dt = A T + B dT/dx + C d/dx(D dT/dx) + E
!!
!! We assume the user is providing the type of boundary conditions and that
!! any Neumann boundary conditions are interpreted as diffs in state variable
!! to be solved (this avoids needing solvers to have to do additional calculations
!! using auxiliary variables not in scope.  

real(wp), dimension(:), intent(in) :: A,B,C,D,E
real(wp), dimension(:), intent(in) :: Ts
real(wp), intent(in) :: Tsminx1, Tsmaxx1, dt
integer, dimension(2), intent(in) :: BCtype  !=0 dirichlet; =1 neumann
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference
integer, parameter :: ll=2                   !number of lower diagonals

real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
real(wp), dimension(size(Ts)) :: Dh
integer :: ix1,lx1

real(wp), dimension(size(Ts)) :: TR

real(wp), dimension(size(Ts)) :: TRBDF21D

!> ORGANIZE SIZES AND THERMAL CONDUCTIVITY
lx1=size(Ts)
Dh(1)=0
Dh(2:lx1)=0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
!    TR(:)=Ts(:)/dt+E(:)
!! boundaries to be overwritten later...  This is now done for each grid point in a separate statement


!! ## TR HALF STEP:  DEFINE A MATRIX USING BANDED STORAGE


! ZZZ - check whether diriclet or neumann...
!> MINX1 BOUNDARY (DIRICHLET)
ix1=1
if (BCtype(1)==0) then
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=0
  M(ll+1,ix1+2)=0
  TR(ix1)=Tsminx1
else
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=-1
  M(ll+1,ix1+2)=0
  TR(ix1)=Tsminx1 !used to be 0, but we change it since we are doing the difference. 
end if


!> FIRST INTERIOR GRID POINT
ix1=2

!> ix1-1
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2 &
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2

!> ix1
M(ll+3,ix1)=1/(dt/2)-A(ix1)/2 &
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2

!> ix1+1, super-diag.
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
M(ll+1,ix1+2)=0
TR(ix1)=Ts(ix1)/(dt/2)+E(ix1) &
  -M(ll+4,ix1-1)*Ts(ix1-1) &
  -(-A(ix1)/2+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2)*Ts(ix1) &
  -M(ll+2,ix1+1)*Ts(ix1+1) &
  -M(ll+1,ix1+2)*Ts(ix1+2)


!> INTERIOR GRID POINTS
do concurrent (ix1=3:lx1-2)
  !! do concurrent OK because only indexing already defined things
  M(ll+5,ix1-2) = 0
  !! ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2 &        !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
  M(ll+3,ix1)=1/(dt/2)-A(ix1)/2 &                           !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2
  M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &    !ix1+1, super-diag.
           -B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
  M(ll+1,ix1+2) = 0
  !! ix1+2 grid point
  TR(ix1)=Ts(ix1)/(dt/2)+E(ix1) &
    -M(ll+5,ix1-2)*Ts(ix1-2) &
    -M(ll+4,ix1-1)*Ts(ix1-1) &
    -(-A(ix1)/2+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2)*Ts(ix1) &
    -M(ll+2,ix1+1)*Ts(ix1+1) &
    -M(ll+1,ix1+2)*Ts(ix1+2)
end do


!> LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2)=0
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2 &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
M(ll+3,ix1)=1/(dt/2)-A(ix1)/2 &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
TR(ix1)=Ts(ix1)/(dt/2)+E(ix1) &
  -M(ll+5,ix1-2)*Ts(ix1-2) &
  -M(ll+4,ix1-1)*Ts(ix1-1) &
  -(-A(ix1)/2+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2)*Ts(ix1) &
  -M(ll+2,ix1+1)*Ts(ix1+1)


! ZZZ - check whether dirichlet or neumann...
!> MAXX1 BOUNDARY
ix1=lx1
if (BCtype(2)==0) then
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=0
  M(ll+3,ix1)=1
  TR(ix1)=Tsmaxx1
else
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=-1
  M(ll+3,ix1)=1
  TR(ix1)=Tsmaxx1
end if


!! ### TR HALF STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER

!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,TR,kl=2)



!! ## BDF2 STEP:  DEFINE A MATRIX USING BANDED STORAGE

!ZZZ - check whether D or N
!> MINX1 BOUNDARY (DIRICHLET)
ix1=1
if (BCtype(1)==0) then
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=0
  M(ll+1,ix1+2)=0
  TRBDF21D(ix1)=Tsminx1
else
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=-1
  M(ll+1,ix1+2)=0
  TRBDF21D(ix1)=Tsminx1
end if


!> FIRST INTERIOR GRID POINT
ix1=2
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/(dt/3)-A(ix1) &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+1,ix1+2)=0
TRBDF21D(ix1)=E(ix1) &
  -1/3._wp*Ts(ix1)/(dt/3) &
  +4/3._wp*TR(ix1)/(dt/3)


!> INTERIOR GRID POINTS
do concurrent (ix1=3:lx1-2)
  M(ll+5,ix1-2)=0                                               !ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &        !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+3,ix1)=1/(dt/3)-A(ix1) &                           !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
  M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &    !ix1+1, super-diag.
           -B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+1,ix1+2)=0                                               !ix1+2 grid point
  TRBDF21D(ix1)=E(ix1) &
    -1/3._wp*Ts(ix1)/(dt/3) &
    +4/3._wp*TR(ix1)/(dt/3)
end do


!LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2)=0
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/(dt/3)-A(ix1) &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))
TRBDF21D(ix1)=E(ix1) &
  -1/3._wp*Ts(ix1)/(dt/3) &
  +4/3._wp*TR(ix1)/(dt/3)

!check whether D or N
!> MAXX1 BOUNDARY
ix1=lx1
if (BCtype(2)==0) then
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=0
  M(ll+3,ix1)=1
  TRBDF21D(ix1)=Tsmaxx1
else
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=-1
  M(ll+3,ix1)=1
  TRBDF21D(ix1)=Tsmaxx1
end if


!! ## BDF2 STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER

!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,TRBDF21D,kl=2)

end function TRBDF21D


function backEuler1D(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,BCtype,dx1,dx1i,coeffs,rhs)

!------------------------------------------------------------
!-------SOLVE A 1D DIFFUSION PROBLEM.  IT IS EXPECTED THAT
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  FORM OF EQUATION
!-------SOLVED IS:
!------- dT/dt + A T + B dT/dx + C d/dx(D dT/dx) = E
!------------------------------------------------------------

real(wp), dimension(:), intent(in) :: A,B,C,D,E
real(wp), dimension(:), intent(in) :: Ts
real(wp), intent(in) :: Tsminx1, Tsmaxx1, dt
integer, dimension(2), intent(in) :: BCtype
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference
real(wp), dimension(:,:), intent(inout), optional :: coeffs
!! intent(out)
real(wp), dimension(:), intent(inout), optional :: rhs
!! intent(out)

integer, parameter :: ll=2                   !number of lower diagonals
real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
real(wp), dimension(size(Ts)) :: Dh
real(wp), dimension(size(Ts)) :: backEuler1D
integer :: ix1,lx1

!------------------------------------------------------------
!-------DEFINE A MATRIX USING BANDED STORAGE
!------------------------------------------------------------
lx1=size(Ts)
Dh(1)=0
Dh(2:lx1)=0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
backEuler1D(:)=Ts(:)/dt+E(:)                !boundaries to be overwritten later...


! check whether D or N
!> MINX1 BOUNDARY, Dirichlet BCS
ix1=1
if (BCtype(1)==0) then
  M(ll+3,ix1)=1       !main diagonal denoted temperature at this grid point... 1*Ts,i=Tsminx1
  M(ll+2,ix1+1)=0     !1st super diagonal
  M(ll+1,ix1+2)=0     !2nd super diagonal
  backEuler1D(ix1)=Tsminx1
else
  M(ll+3,ix1)=1       !main diagonal denoted temperature at this grid point... 1*Ts,i=Tsminx1
  M(ll+2,ix1+1)=-1    !1st super diagonal
  M(ll+1,ix1+2)=0     !2nd super diagonal
  backEuler1D(ix1)=0
end if

!> FIRST INTERIOR GRID POINT
ix1=2
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1, sub-diaginal
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/dt-A(ix1) &                                     !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+1,ix1+2)=0


!> INTERIOR GRID POINTS
do concurrent (ix1=3:lx1-2)
  M(ll+5,ix1-2) = 0
  !! ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+3,ix1)=1/dt-A(ix1) &                                     !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
  M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
           -B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+1,ix1+2) = 0
  !! ix1+2 grid point
end do


!> LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2)=0
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/dt-A(ix1) &                                     !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))

! check whether D or N
!> MAXX1 BOUNDARY
ix1=lx1
if (BCtype(2)==0) then
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=0
  M(ll+3,ix1)=1
  backEuler1D(ix1)=Tsmaxx1
else
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=-1
  M(ll+3,ix1)=1
  backEuler1D(ix1)=0
end if

!> in case we want to output the right-hand side of the system; has to be done here before
!   the output argument (which stores rhs) is overwritten by the solution.
if (present(rhs)) then
  rhs(1:lx1)=backEuler1D(1:lx1)
end if

!! ## DO SOME STUFF TO CALL LAPACK'S BANDED SOLVER
!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,backEuler1D,kl=2)

!> this is for if one wants to output the matrix bands for testing purposes
if (present(coeffs)) then
  coeffs(1:lx1,1)=M(ll+2,1:lx1)
  coeffs(1:lx1,2)=M(ll+3,1:lx1)
  coeffs(1:lx1,3)=M(ll+4,1:lx1)
end if

end function backEuler1D

end module PDEparabolic
