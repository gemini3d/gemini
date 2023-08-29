program test_diffusion1d
!! Need program statement for FORD
!! Solve a time-dependent heat equation in 1D.  See GEMINI-docs repo for
!! a description of the specific problem solved here

use phys_consts, only : wp,pi
use PDEparabolic, only : backEuler1D,TRBDF21D
use h5fortran, only : hdf5_file

implicit none (type, external)

type(hdf5_file) :: h5f
integer, parameter :: npts=256,lt=20*5

character(:), allocatable :: outfile
character(2048) :: argv

character(4) :: ic

real(wp), dimension(npts) :: dx1i
real(wp), dimension(-1:npts+2) :: x1,TsEuler,TsBDF2,Tstrue
real(wp), dimension(npts) :: lambda,A,B,C,D,E
real(wp), dimension(npts+1) :: x1i
real(wp), dimension(0:npts+2) :: dx1
integer :: lx1,it,ix1,ierr
real(wp) :: t=0,dt
real(wp) :: Tsminx1,Tsmaxx1

real(wp), dimension(npts) :: errorEuler,errorBDF2
real(wp), dimension(npts,3) :: abc
real(wp), dimension(npts) :: y

!! create a grid for the calculation
x1=[ (real(ix1-1,wp)/real(npts-1,wp), ix1=-1,npts+2) ]
lx1=npts   !exclude ghost cells in count
dx1=x1(0:lx1+2)-x1(-1:lx1+1)
x1i(1:lx1+1)=0.5*(x1(0:lx1)+x1(1:lx1+1))
dx1i=x1i(2:lx1+1)-x1i(1:lx1)

call get_command_argument(1, argv, status=ierr)
if(ierr/=0) error stop 'please specify input filename'
outfile = trim(argv)

!! write the time, space length and spatial grid to a file
! open(newunit=u,file=outfile,status='replace')
! write(u,*) lt
! write(u,*) lx1
! call writearray(u,x1)
call h5f%open(outfile, action='w')
call h5f%write('/lt', lt)
call h5f%write('/lx1', lx1)
call h5f%write('/x1', x1)

!! initial conditions
TsEuler(-1:lx1+2)=sin(2*pi*x1(-1:lx1+2))+sin(8*pi*x1(-1:lx1+2))
TsBDF2=TsEuler
lambda(:)=1     !thermal conductivity


!! typical diffusion time, make our time step a fraction of this
dt=0.05/8**2/pi**2/maxval(lambda)


!! time iterations
do it=1,lt
  write(ic, '(I4.4)') it
  !boundary values
  Tsminx1=0.0
  Tsmaxx1=0.0

  !solve using two different numerical schemes
  A(:)=0
  B(:)=0
  C(:)=1
  D(:)=lambda(:)
  E(:)=0.0
  
  TsEuler(1:lx1)=backEuler1D(TsEuler(1:lx1),A,B,C,D,E,Tsminx1,Tsmaxx1,dt,[0,0],dx1,dx1i,coeffs=abc,rhs=y)
  TsBDF2(1:lx1)=TRBDF21D(TsBDF2(1:lx1),A,B,C,D,E,Tsminx1,Tsmaxx1,dt,[0,0],dx1,dx1i)
  t=t+dt

  !compute analytical solution to compare
  Tstrue(1:lx1) = exp(-4 *pi**2*lambda*t)*sin(2 *pi*x1(1:lx1))+exp(-64 *pi**2*lambda*t)*sin(8 *pi*x1(1:lx1))

  !! output
  ! write(u,*) t
  ! call writearray(u,TsEuler(1:lx1))
  ! call writearray(u,TsBDF2(1:lx1))
  ! call writearray(u,Tstrue(1:lx1))
  call h5f%write('/t'//ic, t)
  call h5f%write('/TsEuler'//ic, TsEuler(1:lx1))
  call h5f%write('/TsBDF2'//ic, TsBDF2(1:lx1))
  call h5f%write('/TsTrue'//ic, Tstrue(1:lx1))
  call h5f%write('/coeffs'//ic,abc)
  call h5f%write('/rhs'//ic,y)

  !check the validity of the numerical solutions at this time step
  errorEuler(1:lx1)=TsEuler(1:lx1)-Tstrue(1:lx1)
  errorBDF2(1:lx1)=TsBDF2(1:lx1)-Tstrue(1:lx1)
  if (mod(it,5) == 0) then
    print*, 'At time step:  ',it,' max error:  ',maxval(abs(errorEuler)),maxval(abs(errorBDF2))
  end if
  if (maxval(abs(errorEuler)) > 0.05_wp) then    !more that 5% error at a point means something really bad happened...
    print*, 'Time step:  ',it,dt
    error stop 'Excessive error (large max diff) in backward Euler solution, check time step maybe???'
  end if
  if (maxval(abs(errorBDF2)) > 0.05_wp) then
    print*, 'Time step:  ',it,dt
    error stop 'Excessive error (large max diff) in TRBDF2 solution, check time step maybe???'
  end if
end do

! close(u)
call h5f%close()

contains

! subroutine writearray(u,array)
! integer, intent(in) :: u
! real(wp), dimension(:), intent(in) :: array

! integer :: k

! do k=1,size(array)
!   write(u,*) array(k)
! end do
! end subroutine writearray

!
!  subroutine write2Darray(u,array)
!    integer, intent(in) :: u
!    real(wp), dimension(:,:), intent(in) :: array
!
!    integer :: k1,k2
!
!    do k1=1,size(array,1)
!      write(u,'(f8.0)') (array(k1,k2), k2=1,size(array,2))
!    end do
!  end subroutine write2Darray
!

end program
