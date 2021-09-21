module neutral

use phys_consts, only: wp, lnchem, pi, Re, debug
use grid, only: lx1, lx2, lx3
use meshobj, only: curvmesh
use timeutils, only : find_lastdate
use mpimod, only: mpi_cfg
use config, only: gemini_cfg
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use neutraldata2Daxisymmobj, only: neutraldata2Daxisymm
use neutraldata2Dcartobj, only: neutraldata2Dcart

! also links MSIS from vendor/msis00/

implicit none (type, external)
private
public :: Tnmsis, neutral_atmos, make_dneu, clear_dneu, neutral_perturb, neutral_update, init_neutrals, &
  neutral_winds, rotate_geo2native, neutral_denstemp_update, neutral_wind_update, store_geo2native_projections


interface !< atmos.f90
  module subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,nn,Tn,msis_version)
    integer, intent(in) :: ymd(3), msis_version
    real(wp), intent(in) :: UTsecd
    real(wp), dimension(:,:,:), intent(in) :: glat,glon,alt
    real(wp), intent(in) :: activ(3)
    real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3),lnchem), intent(inout) :: nn
    !! intent(out)
    real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)), intent(inout) :: Tn
    !! intent(out)
  end subroutine neutral_atmos
end interface

interface !< perturb.f90
  module subroutine neutral_perturb(cfg,dt,dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dt,dtneu
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    !! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec

    class(curvmesh), intent(inout) :: x
    !! grid structure  (inout because we want to be able to deallocate unit vectors once we are done with them)
    real(wp), intent(in) :: v2grid,v3grid
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    !! intent(out)
    !! neutral params interpolated to plasma grid at requested time
    real(wp), dimension(:,:,:), intent(inout) :: Tn,vn1,vn2,vn3
    !! intent(out)
  end subroutine neutral_perturb
end interface

interface !< wind.f90
  module subroutine neutral_winds(ymd, UTsec, Ap, x, v2grid, v3grid, vn1, vn2, vn3)
    integer, intent(in) :: ymd(3)
    real(wp), intent(in) :: UTsec, Ap
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: v2grid,v3grid
    real(wp), dimension(1:size(x%alt,1),1:size(x%alt,2),1:size(x%alt,3)), intent(inout) :: vn1,vn2,vn3
  end subroutine neutral_winds
end interface

! flag to check whether to apply neutral perturbations
logical :: flagneuperturb=.false.

!! ALL ARRAYS THAT FOLLOW ARE USED WHEN INCLUDING NEUTRAL PERTURBATIONS FROM ANOTHER MODEL
!! ARRAYS TO STORE THE NEUTRAL GRID INFORMATION
!! as long as the neutral module is in scope these persist and do not require a "save"; this variable only used by the axisymmetric interpolation
real(wp), dimension(:), allocatable, target, private :: rhon     !used for axisymmetric 2D simulations, aliased by pointers
real(wp), dimension(:), allocatable, target, private :: yn    !used in cartesian 2D and 3D interpolation
real(wp), dimension(:), allocatable, private :: zn
real(wp), dimension(:), allocatable, private :: xn    !for 3D cartesian interpolation
integer, private :: lrhon,lzn,lyn,lxn

!! STORAGE FOR NEUTRAL SIMULATION DATA.
! These will be singleton in the second dimension (longitude) in the case of 2D interpolation...
!! THESE ARE INCLUDED AS MODULE VARIATIONS TO AVOID HAVING TO REALLOCATE AND DEALLOCIATE EACH TIME WE NEED TO INTERP
real(wp), dimension(:,:,:), allocatable, private :: dnO,dnN2,dnO2,dvnrho,dvnz,dvnx,dTn

!!full grid parameters for root to store input from files.
real(wp), dimension(:), allocatable, private :: xnall
real(wp), dimension(:), allocatable, private :: ynall
integer, private :: lxnall,lynall

!! FIXME: really do not need to store these, particularly all at once...
real(wp), dimension(:,:,:), allocatable, private :: dnOall,dnN2all,dnO2all,dvnrhoall,dvnzall,dvnxall,dTnall

!ARRAYS TO STORE NEUTRAL DATA THAT HAS BEEN INTERPOLATED
real(wp), dimension(:,:,:), allocatable, private :: dnOiprev,dnN2iprev,dnO2iprev,dvnrhoiprev,dvnziprev,dTniprev, &
                                                   dvn1iprev,dvn2iprev,dvn3iprev,dvnxiprev
real(wp), private :: tprev
integer, dimension(3), private :: ymdprev
!! time corresponding to "prev" interpolated data

real(wp), private :: UTsecprev
real(wp), dimension(:,:,:), allocatable, private :: dnOinext,dnN2inext,dnO2inext,dvnrhoinext,dvnzinext, &
                                                   dTninext,dvn1inext,dvn2inext,dvn3inext,dvnxinext
real(wp), private :: tnext
integer, dimension(3), private :: ymdnext
real(wp), private :: UTsecnext

!! data at current time level, (centered in time between current time step and next)
real(wp), dimension(:,:,:), allocatable, protected :: dnOinow,dnN2inow,dnO2inow,dTninow,dvn1inow,dvn2inow,dvn3inow

!SPACE TO STORE PROJECTION FACTORS (rotate from magnetic UEN to curv. dipole...)
real(wp), dimension(:,:,:), allocatable, private :: proj_erhop_e1,proj_ezp_e1,proj_erhop_e2,proj_ezp_e2,proj_erhop_e3,proj_ezp_e3    !these projections are used in the axisymmetric interpolation
real(wp), dimension(:,:,:), allocatable, private :: proj_eyp_e1,proj_eyp_e2,proj_eyp_e3    !these are for Cartesian projections
real(wp), dimension(:,:,:), allocatable, private :: proj_exp_e1,proj_exp_e2,proj_exp_e3

!PLASMA GRID ZI AND RHOI LOCATIONS FOR INTERPOLATIONS
real(wp), dimension(:), allocatable, private :: zi,xi    !this is to be a flat listing of sites on the, rhoi only used in axisymmetric and yi only in cartesian
real(wp), dimension(:), allocatable, target, private :: yi,rhoi

!USED FOR 3D INTERPOLATION WHERE WORKER DIVISIONS ARE COMPLICATED (note that the first dim starts at zero so it matches mpi ID)
real(wp), dimension(:,:), private, allocatable :: extents    !roots array that is used to store min/max x,y,z of each works
integer, dimension(:,:), private, allocatable :: indx       !roots array that contain indices for each workers needed piece of the neutral data
integer, dimension(:,:), private, allocatable :: slabsizes

!! BASE MSIS ATMOSPHERIC STATE ON WHICH TO APPLY PERTURBATIONS
real(wp), dimension(:,:,:,:), allocatable, protected :: nnmsis
real(wp), dimension(:,:,:), allocatable, protected :: Tnmsis
real(wp), dimension(:,:,:), allocatable, protected :: vn1base,vn2base,vn3base

!! projection factors for converting vectors mag->geo; e.g. defining rotation matrix from geographic coords into
real(wp), dimension(:,:,:), allocatable :: proj_ealt_e1,proj_ealt_e2,proj_ealt_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglat_e1,proj_eglat_e2,proj_eglat_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglon_e1,proj_eglon_e2,proj_eglon_e3

!! new module variables for OO refactor
class(neutraldata), allocatable :: atmosperturb

contains

!> initializes neutral atmosphere by:
!    1)  allocating storage space
!    2)  establishing initial background for density, temperature, and winds
!    3)  priming file input so that we have an initial perturbed state to start from (necessary for restart)
subroutine init_neutrals(dt,t,cfg,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
  real(wp), intent(in) :: dt,t
  type(gemini_cfg), intent(in) :: cfg
  integer, dimension(3), intent(in) :: ymd
  real(wp), intent(in) :: UTsec
  class(curvmesh), intent(inout) :: x    ! unit vecs may be deallocated after first setup
  real(wp), intent(in) :: v2grid,v3grid
  real(wp), dimension(:,:,:,:), intent(inout) :: nn
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: Tn
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
  !! intent(out)

  integer, dimension(3) :: ymdtmp
  real(wp) :: UTsectmp
  real(wp) :: tstart,tfin
  
  !! allocation neutral module scope variables so there is space to store all the file input and do interpolations
  call make_dneu()
  
  !! call msis to get an initial neutral background atmosphere
  if (mpi_cfg%myid == 0) call cpu_time(tstart)
  call neutral_atmos(cfg%ymd0,cfg%UTsec0,x%glat,x%glon,x%alt,cfg%activ,nn,Tn,cfg%msis_version)
  if (mpi_cfg%myid == 0) then
    call cpu_time(tfin)
    print *, 'Initial neutral density and temperature (from MSIS) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
  end if
  
  !> Horizontal wind model initialization/background
  if (mpi_cfg%myid == 0) call cpu_time(tstart)
  call neutral_winds(cfg%ymd0, cfg%UTsec0, Ap=cfg%activ(3), x=x, v2grid=v2grid,v3grid=v3grid,vn1=vn1, vn2=vn2, vn3=vn3)
  !! we sum the horizontal wind with the background state vector
  !! if HWM14 is disabled, neutral_winds returns the background state vector unmodified
  if (mpi_cfg%myid == 0) then
    call cpu_time(tfin)
    print *, 'Initial neutral winds (from HWM) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
  end if


  !! perform an initialization for the perturbation quantities
  if (cfg%flagdneu==1) then
    ! set flag denoted neutral perturbations
    flagneuperturb=.true.

    ! allocate correct type, FIXME: eventuallly no shunt to 3D
    select case (cfg%interptype)
    case (0)
      allocate(neutraldata2Dcart::atmosperturb)
    case (1)
      allocate(neutraldata2Daxisymm::atmosperturb)
    case (3)
      allocate(neutraldata3D::atmosperturb)
    case default
      error stop 'non-standard neutral interpolation type chosen in config.nml...'
    end select

    ! call object init procedure
    call atmosperturb%init(cfg,cfg%sourcedir,x,dt,cfg%dtneu,ymd,UTsec)
  end if
end subroutine init_neutrals


!> update density, temperature, and winds
subroutine neutral_update(nn,Tn,vn1,vn2,vn3,v2grid,v3grid)
  !! adds stored base and perturbation neutral atmospheric parameters
  !!  these are module-scope parameters so not needed as input
  real(wp), dimension(:,:,:,:), intent(inout) :: nn
  !! intent(out)  
  real(wp), dimension(:,:,:), intent(inout) :: Tn
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
  !! intent(out)
  real(wp) :: v2grid,v3grid

  call neutral_denstemp_update(nn,Tn)
  call neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid)
end subroutine neutral_update


!> Adds stored base (viz. background) and perturbation neutral atmospheric density
subroutine neutral_denstemp_update(nn,Tn)
  real(wp), dimension(:,:,:,:), intent(out) :: nn
  real(wp), dimension(:,:,:), intent(out) :: Tn
  
  !> background neutral parameters
  nn=nnmsis
  Tn=Tnmsis
  
  !> add perturbations, if used
  if (flagneuperturb) then
    nn(:,:,:,1)=nn(:,:,:,1)+dnOinow
    nn(:,:,:,2)=nn(:,:,:,2)+dnN2inow
    nn(:,:,:,3)=nn(:,:,:,3)+dnO2inow
    nn(:,:,:,1)=max(nn(:,:,:,1),1._wp)
    nn(:,:,:,2)=max(nn(:,:,:,2),1._wp)
    nn(:,:,:,3)=max(nn(:,:,:,3),1._wp)
    !! note we are not adjusting derived densities like NO since it's not clear how they may be related to major
    !! species perturbations.
  
    Tn=Tn+dTninow
    Tn=max(Tn,50._wp)
  end if
end subroutine neutral_denstemp_update


!> update wind variables with background and perturbation quantities
subroutine neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid)
  real(wp), dimension(:,:,:), intent(out) :: vn1,vn2,vn3
  real(wp) :: v2grid,v3grid
  
  !> background neutral parameters
  vn1=vn1base
  vn2=vn2base
  vn3=vn3base
  
  !> perturbations, if used
  if (flagneuperturb) then
    vn1=vn1+dvn1inow
    vn2=vn2+dvn2inow
    vn3=vn3+dvn3inow
  end if
  
  !> subtract off grid drift speed (needs to be set to zero if not lagrangian grid)
  vn2=vn2-v2grid
  vn3=vn3-v3grid
end subroutine neutral_wind_update


!> rotate winds from geographic to model native coordinate system (x1,x2,x3)
subroutine rotate_geo2native(vnalt,vnglat,vnglon,x,vn1,vn2,vn3)
  real(wp), dimension(:,:,:), intent(in) :: vnalt,vnglat,vnglon
  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:), intent(out) :: vn1,vn2,vn3
  real(wp), dimension(1:size(vnalt,1),1:size(vnalt,2),1:size(vnalt,3),3) :: ealt,eglat,eglon
  integer :: lx1,lx2,lx3

  !> if first time called then allocate space for projections and compute
  if (.not. allocated(proj_ealt_e1)) then
    call x%calc_unitvec_geo(ealt,eglon,eglat)
    call store_geo2native_projections(x,ealt,eglon,eglat)
  end if

  !> rotate vectors into model native coordinate system
  vn1=vnalt*proj_ealt_e1+vnglat*proj_eglat_e1+vnglon*proj_eglon_e1
  vn2=vnalt*proj_ealt_e2+vnglat*proj_eglat_e2+vnglon*proj_eglon_e2
  vn3=vnalt*proj_ealt_e3+vnglat*proj_eglat_e3+vnglon*proj_eglon_e3
end subroutine rotate_geo2native


!> compute projections for rotating winds geographic to native coordinate system
subroutine store_geo2native_projections(x,ealt,eglon,eglat,rotmat)
  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:,:), intent(in) :: ealt,eglon,eglat
  real(wp), dimension(:,:,:,:,:), intent(out), optional :: rotmat    ! for debugging purposes
  integer :: ix1,ix2,ix3,lx1,lx2,lx3

  !! allocate module-scope space for the projection factors
  lx1=size(ealt,1); lx2=size(ealt,2); lx3=size(ealt,3);
  allocate(proj_ealt_e1(lx1,lx2,lx3),proj_eglat_e1(lx1,lx2,lx3),proj_eglon_e1(lx1,lx2,lx3))
  allocate(proj_ealt_e2(lx1,lx2,lx3),proj_eglat_e2(lx1,lx2,lx3),proj_eglon_e2(lx1,lx2,lx3))
  allocate(proj_ealt_e3(lx1,lx2,lx3),proj_eglat_e3(lx1,lx2,lx3),proj_eglon_e3(lx1,lx2,lx3))

  !! compute projections (dot products of unit vectors)
  proj_ealt_e1=sum(ealt*x%e1,4)
  proj_eglat_e1=sum(eglat*x%e1,4)
  proj_eglon_e1=sum(eglon*x%e1,4)
  proj_ealt_e2=sum(ealt*x%e2,4)
  proj_eglat_e2=sum(eglat*x%e2,4)
  proj_eglon_e2=sum(eglon*x%e2,4)
  proj_ealt_e3=sum(ealt*x%e3,4)
  proj_eglat_e3=sum(eglat*x%e3,4)
  proj_eglon_e3=sum(eglon*x%e3,4)

  !! store the rotation matrix to convert geo to native if the user wants it
  if (present(rotmat)) then
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          rotmat(1,1:3,ix1,ix2,ix3)=[proj_ealt_e1(ix1,ix2,ix3),proj_eglat_e1(ix1,ix2,ix3),proj_eglon_e1(ix1,ix2,ix3)]
          rotmat(2,1:3,ix1,ix2,ix3)=[proj_ealt_e2(ix1,ix2,ix3),proj_eglat_e2(ix1,ix2,ix3),proj_eglon_e2(ix1,ix2,ix3)]
          rotmat(3,1:3,ix1,ix2,ix3)=[proj_ealt_e3(ix1,ix2,ix3),proj_eglat_e3(ix1,ix2,ix3),proj_eglon_e3(ix1,ix2,ix3)]
        end do
      end do
    end do
  end if
end subroutine store_geo2native_projections


subroutine make_dneu()
!ZZZ - could make this take in type of neutral interpolation and do allocations accordingly

!allocate and compute plasma grid z,rho locations and space to save neutral perturbation variables and projection factors
allocate(zi(lx1*lx2*lx3),rhoi(lx1*lx2*lx3))
allocate(yi(lx1*lx2*lx3))
allocate(xi(lx1*lx2*lx3))
allocate(proj_erhop_e1(lx1,lx2,lx3),proj_ezp_e1(lx1,lx2,lx3),proj_erhop_e2(lx1,lx2,lx3),proj_ezp_e2(lx1,lx2,lx3), &
         proj_erhop_e3(lx1,lx2,lx3),proj_ezp_e3(lx1,lx2,lx3))
allocate(proj_eyp_e1(lx1,lx2,lx3),proj_eyp_e2(lx1,lx2,lx3),proj_eyp_e3(lx1,lx2,lx3))
allocate(proj_exp_e1(lx1,lx2,lx3),proj_exp_e2(lx1,lx2,lx3),proj_exp_e3(lx1,lx2,lx3))
allocate(dnOiprev(lx1,lx2,lx3),dnN2iprev(lx1,lx2,lx3),dnO2iprev(lx1,lx2,lx3),dvnrhoiprev(lx1,lx2,lx3), &
         dvnziprev(lx1,lx2,lx3),dTniprev(lx1,lx2,lx3),dvn1iprev(lx1,lx2,lx3),dvn2iprev(lx1,lx2,lx3), &
         dvn3iprev(lx1,lx2,lx3))
allocate(dvnxiprev(lx1,lx2,lx3))
allocate(dnOinext(lx1,lx2,lx3),dnN2inext(lx1,lx2,lx3),dnO2inext(lx1,lx2,lx3),dvnrhoinext(lx1,lx2,lx3), &
         dvnzinext(lx1,lx2,lx3),dTninext(lx1,lx2,lx3),dvn1inext(lx1,lx2,lx3),dvn2inext(lx1,lx2,lx3), &
         dvn3inext(lx1,lx2,lx3))
allocate(dvnxinext(lx1,lx2,lx3))
allocate(nnmsis(lx1,lx2,lx3,lnchem),Tnmsis(lx1,lx2,lx3),vn1base(lx1,lx2,lx3),vn2base(lx1,lx2,lx3),vn3base(lx1,lx2,lx3))
allocate(dnOinow(lx1,lx2,lx3),dnN2inow(lx1,lx2,lx3),dnO2inow(lx1,lx2,lx3),dvn1inow(lx1,lx2,lx3),dvn2inow(lx1,lx2,lx3), &
           dvn3inow(lx1,lx2,lx3), dTninow(lx1,lx2,lx3))

!start everyone out at zero
zi = 0
rhoi = 0
yi = 0
xi = 0
proj_erhop_e1 = 0
proj_ezp_e1 = 0
proj_erhop_e2 = 0
proj_ezp_e2 = 0
proj_erhop_e3 = 0
proj_ezp_e3 = 0
proj_eyp_e1 = 0
proj_eyp_e2 = 0
proj_eyp_e3 = 0
proj_exp_e1 = 0
proj_exp_e2 = 0
proj_exp_e3 = 0
dnOiprev = 0
dnN2iprev = 0
dnO2iprev = 0
dTniprev = 0
dvnrhoiprev = 0
dvnziprev = 0
dvn1iprev = 0
dvn2iprev = 0
dvn3iprev = 0
dvnxiprev = 0
dnOinext = 0
dnN2inext = 0
dnO2inext = 0
dTninext = 0
dvnrhoinext = 0
dvnzinext = 0
dvn1inext = 0
dvn2inext = 0
dvn3inext = 0
dvnxinext = 0
nnmsis = 0
Tnmsis = 0
vn1base = 0
vn2base = 0
vn3base = 0
dnOinow = 0
dnN2inow = 0
dnO2inow = 0
dTninow = 0
dvn1inow = 0
dvn2inow = 0
dvn3inow = 0

!now initialize some module variables
tprev = 0
tnext = 0

end subroutine make_dneu


subroutine clear_dneu

!stuff allocated at beginning of program
deallocate(zi,rhoi)
deallocate(yi)
deallocate(proj_erhop_e1,proj_ezp_e1,proj_erhop_e2,proj_ezp_e2, &
         proj_erhop_e3,proj_ezp_e3)
deallocate(proj_eyp_e1,proj_eyp_e2,proj_eyp_e3)
deallocate(dnOiprev,dnN2iprev,dnO2iprev,dvnrhoiprev,dvnziprev,dTniprev,dvn1iprev,dvn2iprev,dvn3iprev)
deallocate(dnOinext,dnN2inext,dnO2inext,dvnrhoinext,dvnzinext,dTninext,dvn1inext,dvn2inext,dvn3inext)
deallocate(nnmsis,Tnmsis,vn1base,vn2base,vn3base)
deallocate(dnOinow,dnN2inow,dnO2inow,dTninow,dvn1inow,dvn2inow,dvn3inow)

!check whether any other module variables were allocated and deallocate accordingly
if (allocated(zn) ) then    !if one is allocated, then they all are
  deallocate(zn)
  deallocate(dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
end if
if (allocated(rhon)) then
  deallocate(rhon)
end if
if (allocated(yn)) then
  deallocate(yn)
end if
if (allocated(extents)) then
  deallocate(extents,indx,slabsizes)
end if
if (allocated(dvnx)) then
  deallocate(dvnx)
end if
if (allocated(xn)) then
  deallocate(xn)
end if
if (allocated(xnall)) then
  deallocate(xnall,ynall)
end if
if (allocated(dnOall)) then    !! 3D input doesn't allocate this so check independent of coord-all
  deallocate(dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
end if
if (allocated(proj_ealt_e1)) then
  deallocate(proj_ealt_e1,proj_eglat_e1,proj_eglon_e1)
  deallocate(proj_ealt_e2,proj_eglat_e2,proj_eglon_e2)
  deallocate(proj_ealt_e3,proj_eglat_e3,proj_eglon_e3)
end if

end subroutine clear_dneu

end module neutral
