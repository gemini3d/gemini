program fullgrid_testdriver

use, intrinsic :: ISO_Fortran_env,  only : wp=>real64
use dipole, only : make_dipolemesh,dipolemesh,de_dipolemesh

implicit none

integer, parameter :: lq=256,lp=192,lphi=128
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(lphi) :: phi
real(wp), dimension(2), parameter :: qlims=[-0.5851937,0.5851937]
real(wp), dimension(2), parameter :: plims=[1.2053761,1.5820779]
real(wp), dimension(2), parameter :: philims=[2.0,2.5]
integer :: iq,ip,iphi
type(dipolemesh) :: x
real(wp) :: minchkvar,maxchkvar


! define a grid, in reality this would be pull in from a file
q=[(qlims(1) + (qlims(2)-qlims(1)/lq*(iq-1)),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1)/lp*(ip-1)),ip=1,lp)]
phi=[(philims(1) + (philims(2)-philims(1))/lphi*(iphi-1),iphi=1,lphi)]

! call grid generation for this grid def.
print*, 'fullgrid_testdriver:  Calling dipole mesh constructor...'
x=make_dipolemesh(q,p,phi)

! now do some basic sanity checks
print*, 'fullgrid_testdriver:  Starting basic checks...'
minchkvar=minval(x%er); maxchkvar=maxval(x%er);
print*, ' fullgrid_testdriver, er:  ',minchkvar,maxchkvar
minchkvar=minval(x%etheta); maxchkvar=maxval(x%ephi);
print*, ' fullgrid_testdriver, etheta:  ',minchkvar,maxchkvar
minchkvar=minval(x%ephi); maxchkvar=maxval(x%ephi);
print*, ' fullgrid_testdriver, ephi:  ',minchkvar,maxchkvar
minchkvar=minval(x%eq); maxchkvar=maxval(x%eq);
print*, ' fullgrid_testdriver, eq:  ',minchkvar,maxchkvar
minchkvar=minval(x%ep); maxchkvar=maxval(x%ep);
print*, ' fullgrid_testdriver, ep:  ',minchkvar,maxchkvar
minchkvar=minval(x%Bmag); maxchkvar=maxval(x%Bmag);
print*, ' fullgrid_testdriver, Bmag (nT):  ',minchkvar*1e9,maxchkvar*1e9

! deallocate the grid before ending the program
print*, 'fullgrid_testdriver:  Deallocating mesh...'
call de_dipolemesh(x)

end program fullgrid_testdriver