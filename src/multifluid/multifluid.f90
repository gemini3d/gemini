! Copyright 2021 Matthew Zettergren

! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!>  This module contains subroutines to be called by a gemini main program in order to execute numerical solutions
!     to the ionospheric transport equations (electrodynamics calls are located elsewhere).
module multifluid

use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
use advec, only: interface_vels_allspec,sweep3_allspec,sweep1_allspec,sweep2_allspec
use calculus, only: etd_uncoupled, div3d
use collisions, only:  thermal_conduct
use phys_consts, only : wp,pi,qs,lsp,gammas,kB,ms,mindensdiv,mindens,mindensnull, debug
use diffusion, only:  trbdf23d, diffusion_prep, backEuler3D
use grid, only: lx1, lx2, lx3, gridflag
use meshobj, only: curvmesh
use ionization, only: ionrate_glow98, ionrate_fang, eheating, photoionization
use precipBCs_mod, only: precipBCs_fileinput, precipBCs
use sources, only: srcsenergy, srcsmomentum, srcscontinuity
use timeutils, only : sza
use gemini3d_config, only: gemini_cfg
use precipdataobj, only: precipdata

implicit none (type, external)
private
public :: sweep3_allparams,sweep1_allparams,sweep2_allparams,source_loss_allparams,VNRicht_artvisc,compression, &
            energy_diffusion,impact_ionization,clean_param,rhoe2T,T2rhoe,rhov12v1,v12rhov1,clean_param_after_regrid

integer, parameter :: lprec=2
!! number of precipitating electron populations
real(wp), allocatable, dimension(:,:,:,:) :: PrPrecipG
real(wp), allocatable, dimension(:,:,:) :: QePrecipG, iverG

real(wp), parameter :: xicon = 3
!real(wp), parameter :: xicon = 0
!! artificial viscosity, decent value for closed field-line grids extending to high altitudes, can be set to 0 for cartesian simulations not exceed altitudes of 1500 km.

contains
!> sweep advection for all plasma parameters in the x3 direction
subroutine sweep3_allparams(dt,x,vs3i,ns,rhovs1,rhoes)
  real(wp), intent(in) :: dt
  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:,:), intent(in) :: vs3i
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,rhoes
!  if (minval(rhoes) < 0) then
!    print*, '1 rhoes data trashed:  ',minval(rhoes),maxval(rhoes),minloc(rhoes),maxloc(rhoes)
!    print*, 'vs3i:  ',minval(vs3i),maxval(vs3i),minloc(vs3i),maxloc(vs3i)
!    print*, 'ns:  ',minval(ns),maxval(ns),minloc(ns),maxloc(ns)
!    open(newunit=funit,file='error.dat',status='replace',access='stream')
!    write(funit) ns
!    write(funit) vs3i
!    close(funit)
!    error stop
!  end if

  call sweep3_allspec(ns,vs3i,dt,x,0,6)
  call sweep3_allspec(rhovs1,vs3i,dt,x,1,6)
  call sweep3_allspec(rhoes,vs3i,dt,x,0,7)

!  if (minval(rhoes) < 0) then
!    print*, '2 rhoes data trashed:  ',minval(rhoes),maxval(rhoes),minloc(rhoes),maxloc(rhoes)
!    print*, 'vs3i:  ',minval(vs3i),maxval(vs3i),minloc(vs3i),maxloc(vs3i)
!    print*, 'ns:  ',minval(ns),maxval(ns),minloc(ns),maxloc(ns)
!    open(newunit=funit,file='error.dat',status='replace',access='stream')
!    write(funit) ns
!    write(funit) vs3i
!    close(funit)
!    error stop
!  end if
end subroutine sweep3_allparams


!> sweep all parameters in the x1 direction
subroutine sweep1_allparams(dt,x,vs1i,ns,rhovs1,rhoes)
  real(wp), intent(in) :: dt
  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:,:), intent(in) :: vs1i
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,rhoes

  call sweep1_allspec(ns,vs1i,dt,x,6)     ! sweep1 doesn't need to know the rank of the advected quantity
  call sweep1_allspec(rhovs1,vs1i,dt,x,6)
  call sweep1_allspec(rhoes,vs1i,dt,x,7)
end subroutine sweep1_allparams


!> sweep all parameters in the x2 direction
subroutine sweep2_allparams(dt,x,vs2i,ns,rhovs1,rhoes)
  real(wp), intent(in) :: dt
  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:,:), intent(in) :: vs2i
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,rhoes

  call sweep2_allspec(ns,vs2i,dt,x,0,6)
  call sweep2_allspec(rhovs1,vs2i,dt,x,1,6)
  call sweep2_allspec(rhoes,vs2i,dt,x,0,7)
end subroutine sweep2_allparams


!> execute diffusion of energy and then source/loss terms for all equations
subroutine source_loss_allparams(dt,t,cfg,ymd,UTsec,x,E1,E2,E3,Q,f107a,f107,nn,vn1,vn2,vn3, &
                                   Tn,first,ns,rhovs1,rhoes,vs1,vs2,vs3,Ts,iver,gavg,Tninf, &
                                   eprecip)
  real(wp), intent(in) :: dt,t
  type(gemini_cfg), intent(in) :: cfg
  integer, dimension(3), intent(in) :: ymd
  real(wp), intent(in) :: UTsec
  class(curvmesh), intent(in) :: x
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: E1, E2, E3
  real(wp), dimension(:,:,:,:), intent(in) :: Q
  real(wp), intent(in) :: f107a,f107
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
  logical, intent(in) :: first
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,rhoes,vs1,vs2,vs3,Ts
  real(wp), dimension(:,:,:), intent(inout) :: iver
  real(wp), intent(in) :: gavg,Tninf
  type(precipdata), intent(inout) :: eprecip
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)) :: Pr,Lo
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)-1) :: Prprecip
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: Qeprecip
  real(wp), dimension(1:size(ns,2)-4,1:size(ns,3)-4,lprec) :: W0,PhiWmWm2
  real(wp) :: tstart,tfin

  !print*, 'Begin src/loss:  ',minval(E1),maxval(E1)


  !> Establish top boundary conditions for electron precipitation
  if (cfg%flagprecfile==1) then
    call precipBCs_fileinput(dt,t,cfg,ymd,UTsec,x,W0,PhiWmWm2,eprecip)
  else
    !! no file input specified, so just call 'regular' function
    call precipBCs(cfg,W0,PhiWmWm2)
  end if

  ! Stiff/balanced energy source, i.e. source/losses for energy equation(s)
  call cpu_time(tstart)
  Prprecip=0.0    ! procedures accumulate rates so need to initialize to zero each time step before rates are updated
  Qeprecip=0.0
  call impact_ionization(cfg,t,dt,x,ymd,UTsec,f107a,f107,Prprecip,Qeprecip,W0,PhiWmWm2,iver,ns,Ts,nn,Tn,first)   ! precipiting electrons
  call solar_ionization(t,x,ymd,UTsec,f107a,f107,Prprecip,Qeprecip,ns,nn,Tn,gavg,Tninf)     ! solar ionization source
  call srcsEnergy(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo,E2,E3,x,cfg)                     ! collisional interactions
  call energy_source_loss(dt,Pr,Lo,Qeprecip,rhoes,Ts,ns)                         ! source/loss numerical solution
  call cpu_time(tfin)
  !if (mpi_cfg%myid==0 .and. debug) then
  !  print *, 'Energy sources substep for time step:  ',t,'done in cpu_time of:  ',tfin-tstart
  !end if

  !print*, 'After energy substep:  ',minval(E1),maxval(E1)

  !ALL VELOCITY SOURCES
  call cpu_time(tstart)
  call srcsMomentum(nn,vn1,Tn,ns,vs1,vs2,vs3,Ts,E1,Q,x,Pr,Lo)    !added artificial viscosity...
  call momentum_source_loss(dt,x,Pr,Lo,ns,rhovs1,vs1)
  call cpu_time(tfin)
  !if (mpi_cfg%myid==0 .and. debug) then
  !  print *, 'Velocity sources substep for time step:  ',t,'done in cpu_time of:  ',tfin-tstart
  !end if

  !ALL MASS SOURCES
  call cpu_time(tstart)
  call srcsContinuity(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo)
  call mass_source_loss(dt,Pr,Lo,Prprecip,ns)
  call cpu_time(tfin)
  !if (mpi_cfg%myid==0 .and. debug) then
  !  print *, 'Mass sources substep for time step:  ',t,'done in cpu_time of:  ',tfin-tstart
  !end if
end subroutine source_loss_allparams


!> Compute temperatures from internal energy densities
subroutine rhoe2T(ns,rhoes,Ts)
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  ns,rhoes,Ts
  integer :: isp,lsp

  lsp=size(rhoes,4)
  do isp=1,lsp
    Ts(:,:,:,isp)=(gammas(isp) - 1)/kB*rhoes(:,:,:,isp)/max(ns(:,:,:,isp),mindensdiv)
    Ts(:,:,:,isp)=max(Ts(:,:,:,isp), 100._wp)
    !! convert internal specific energy density back into temperature
  end do
end subroutine rhoe2T


!> Convert temperature to specific internal energy density
subroutine T2rhoe(ns,Ts,rhoes)
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  ns,rhoes,Ts
  integer :: isp,lsp

  ! convert temperature to enegy density
  lsp=size(Ts,4)
  do isp=1,lsp
    rhoes(:,:,:,isp)=ns(:,:,:,isp)*kB*Ts(:,:,:,isp)/(gammas(isp) - 1)
  end do
end subroutine T2rhoe


!> Convert velocity to momentum
subroutine v12rhov1(ns,vs1,rhovs1)
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  ns,vs1,rhovs1
  integer :: isp,lsp

  lsp=size(vs1,4)
  do isp=1,lsp
    rhovs1(:,:,:,isp)=ns(:,:,:,isp)*ms(isp)*vs1(:,:,:,isp)
  end do
end subroutine v12rhov1


!> Compute electron density and velocity given ion momenta, compute ion velocities as well
subroutine rhov12v1(ns,rhovs1,vs1)
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  ns,rhovs1,vs1
  integer :: isp,lsp
  real(wp), dimension(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2) :: chrgflux

  lsp=size(ns,4)

  chrgflux=0.0
  do isp=1,lsp-1
    vs1(:,:,:,isp)=rhovs1(:,:,:,isp)/(ms(isp)*max(ns(:,:,:,isp),mindensdiv))
    chrgflux=chrgflux+ns(:,:,:,isp)*qs(isp)*vs1(:,:,:,isp)
  end do
  ns(:,:,:,lsp)=sum(ns(:,:,:,1:lsp-1),4)
!!      vs1(1:lx1,1:lx2,1:lx3,lsp)=1/ns(1:lx1,1:lx2,1:lx3,lsp)/qs(lsp)*(J1-chrgflux)   !density floor needed???
  vs1(:,:,:,lsp)=-1/max(ns(:,:,:,lsp),mindensdiv)/qs(lsp)*chrgflux   !really not strictly correct, should include current density
end subroutine rhov12v1


!> Compute artifical viscosity
subroutine VNRicht_artvisc(ns,vs1,Q)
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)), intent(inout) :: Q
  real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4) :: v1iupdate
!! temp interface velocities for art. viscosity
  real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-4) :: dv1iupdate
!! interface diffs. for art. visc.
  integer :: isp,lsp

  lsp=size(ns,4)

  !print*, shape(vs1)
  !print*, shape(ns)

   !print*, 'viscous:  ',shape(vs1(1:lx1,1:lx2,1:lx3,:)),minval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       maxval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       minloc(vs1(1:lx1,1:lx2,1:lx3,:)),maxloc(vs1(1:lx1,1:lx2,1:lx3,:))

!  !ARTIFICIAL VISCOSITY (NOT REALLY NEED BELOW 1000 KM ALT.).  NOTE THAT WE DON'T CHECK WHERE SUBCYCLING IS NEEDED SINCE, IN MY EXPERIENCE THEN CODE IS BOMBING ANYTIME IT IS...
!  ! Interestingly, this is accessing ghost cells of velocity so if they are overwritten by clean_params this viscosity calculation would generate "odd" results
  if (xicon>0) then
    do isp=1,lsp-1
      v1iupdate(1:lx1+1,:,:)=0.5_wp*(vs1(0:lx1,1:lx2,1:lx3,isp)+vs1(1:lx1+1,1:lx2,1:lx3,isp))    !compute an updated interface velocity (only in x1-direction)
      dv1iupdate=v1iupdate(2:lx1+1,:,:)-v1iupdate(1:lx1,:,:)
      Q(:,:,:,isp)=ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*0.25_wp*xicon**2*(min(dv1iupdate,0._wp))**2   !note that viscosity does not have/need ghost cells
    end do
    Q(:,:,:,lsp) = 0
  else
    Q=0._wp
  end if
   !Q=0.0
end subroutine VNRicht_artvisc


!> Adiabatic compression term, including (precomputed) artifical viscosity.  All velocities must be haloed a single
!    point prior to calling this procedure.  Upon entering this procedure the specific internal energy density contains
!    the most recent updated state, while the temperature may or mat not.  Upon exit only the energy density is updated.
subroutine compression(dt,x,vs1,vs2,vs3,Q,rhoes)
  real(wp), intent(in) :: dt
  class(curvmesh), intent(in) :: x
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs1,vs2,vs3
  real(wp), dimension(:,:,:,:), intent(in) :: Q
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: rhoes
  real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-4) :: paramtrim,rhoeshalf
  real(wp), dimension(0:size(vs1,1)-3,0:size(vs1,2)-3,0:size(vs1,3)-3) :: divvs
  integer :: isp,lsp

   !print*, 'compression:  ',  shape(vs1(1:lx1,1:lx2,1:lx3,:)),minval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       maxval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       minloc(vs1(1:lx1,1:lx2,1:lx3,:)),maxloc(vs1(1:lx1,1:lx2,1:lx3,:))

  lsp=size(vs1,4)
  do isp=1,lsp
    divvs = div3D(vs1(0:lx1+1,0:lx2+1,0:lx3+1,isp),&
                  vs2(0:lx1+1,0:lx2+1,0:lx3+1,isp), &
                  vs3(0:lx1+1,0:lx2+1,0:lx3+1,isp),x,0,lx1+1,0,lx2+1,0,lx3+1)
    !! diff with one set of ghost cells to preserve second order accuracy over the grid
    paramtrim=rhoes(1:lx1,1:lx2,1:lx3,isp)

    rhoeshalf = paramtrim - dt/2 * (paramtrim*(gammas(isp)-1) + Q(:,:,:,isp)) * divvs(1:lx1,1:lx2,1:lx3)
    !! t+dt/2 value of internal energy, use only interior points of divvs for second order accuracy

    paramtrim=paramtrim-dt*(rhoeshalf*(gammas(isp) - 1)+Q(:,:,:,isp))*divvs(1:lx1,1:lx2,1:lx3)
    rhoes(1:lx1,1:lx2,1:lx3,isp)=paramtrim
  end do
end subroutine compression


!> Execute energy diffusion substep, no mpi required.  Upon entering this procedure the temperature needs to be have its
!     most recently updated state.  Upon exit the temperature will be updated with diffusion properly applied.
subroutine energy_diffusion(dt,x,ns,Ts,J1,nn,Tn,flagdiffsolve,Teinf)
  real(wp), intent(in) :: dt
  class(curvmesh), intent(in) :: x
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: Ts
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: J1
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: Tn
  integer, intent(in) :: flagdiffsolve
  real(wp), intent(in) :: Teinf
  real(wp), dimension(-1:size(Ts,1)-2,-1:size(Ts,2)-2,-1:size(Ts,3)-2) :: param
  real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4) :: A,B,C,D,E,lambda,beta
  integer :: isp,lsp

  lsp=size(Ts,4)
  do isp=1,lsp
    param=Ts(:,:,:,isp)     !temperature for this species
    call thermal_conduct(isp,param,ns(:,:,:,isp),nn,J1,lambda,beta)

    call diffusion_prep(isp,x,lambda,beta,ns(:,:,:,isp),param,A,B,C,D,E,Tn,Teinf)
    select case (flagdiffsolve)
      case (1)
        param=backEuler3D(param,A,B,C,D,E,dt,x)    !1st order method, only use if you are seeing grid-level oscillations in temperatures
      case (2)
        param=TRBDF23D(param,A,B,C,D,E,dt,x)       !2nd order method, should be used for most simulations
      case default
        print*, 'Unsupported diffusion solver type/mode:  ',flagdiffsolve,'.  Should be either 1 or 2.'
        error stop
    end select

    Ts(:,:,:,isp) = param
    Ts(:,:,:,isp) = max(Ts(:,:,:,isp), 100._wp)    ! is this necessary or does clean_param take care of???
  end do
end subroutine energy_diffusion


!> *Accumulates* ionization and heating rates into Prprecip,Qeprecip arrays; note that if you want only
!    rates from impact ionization these arrays will need to be initialized to zero before calling this
!    procedure.  Note that this procedure does need updated density and temperature data (i.e. ns and Ts)
subroutine impact_ionization(cfg,t,dt,x,ymd,UTsec,f107a,f107,Prprecip,Qeprecip,W0,PhiWmWm2,iver,ns,Ts,nn,Tn,first)
  type(gemini_cfg), intent(in) :: cfg
  real(wp), intent(in) :: t,dt
  class(curvmesh), intent(in) :: x
  integer, dimension(3), intent(in) :: ymd
  real(wp), intent(in) :: UTsec
  real(wp), intent(in) :: f107a,f107
  real(wp), dimension(:,:,:,:), intent(inout) :: Prprecip
  real(wp), dimension(:,:,:), intent(inout) :: Qeprecip
  real(wp), dimension(:,:,:), intent(in) :: W0,PhiWmWm2
  real(wp), dimension(:,:,:), intent(inout) :: iver
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: Tn
  logical, intent(in) :: first  !< first time step
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)-1) :: Prpreciptmp
  integer :: iprec,lprec
  !! note PrprecipG and the like are module-scope variables

  ! first check that our module-scope arrays are allocated before going on to calculations
  if ((cfg%flagglow/=0).and.(.not.allocated(PrprecipG))) then
    allocate(PrprecipG(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)-1))
    PrprecipG(:,:,:,:)=0
  end if
  if ((cfg%flagglow/=0).and.(.not.allocated(QeprecipG))) then
    allocate(QeprecipG(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4))
    QeprecipG(:,:,:)=0
  end if
  if ((cfg%flagglow/=0).and.(.not.allocated(iverG))) then
    allocate(iverG(size(iver,1),size(iver,2),size(iver,3)))
    iverG(:,:,:)=0
  end if

  ! compute impact ionization given input boundary conditions
  lprec=size(W0,3)    ! just recompute the number of precipitating populations
  if (gridflag/=0) then
    if (cfg%flagglow==0) then
      !! Fang et al 2008 parameterization
      do iprec=1,lprec
        !! loop over the different populations of precipitation (2 here?), accumulating production rates
        Prpreciptmp = ionrate_fang(W0(:,:,iprec), PhiWmWm2(:,:,iprec), nn, Tn, cfg%flag_fang, x%g1)
        !! calculation based on Fang et al [2008]
        Prprecip=Prprecip+Prpreciptmp
      end do
      Prprecip = max(Prprecip, 1e-5_wp)         ! should resort to fill values only after all populations accumulated
      Qeprecip = eheating(nn,Prprecip,ns)    ! once we have total ionization rate (all populations) compute the elec. heating rate
    else
      !! glow model
      if (int(t/cfg%dtglow)/=int((t+dt)/cfg%dtglow) .or. first) then
        !if (mpi_cfg%myid==0) print*, 'Note:  preparing to call GLOW...  This could take a while if your grid is large...'
        PrprecipG=0; QeprecipG=0; iverG=0;
        call ionrate_glow98(W0,PhiWmWm2,ymd,UTsec,f107,f107a,x%glat(1,1:lx2,1:lx3),x%glon(1,1:lx2,1:lx3), &
                            x%alt(1:lx1,1:lx2,1:lx3),nn,Tn,ns,Ts, &
                            QeprecipG, iverG, PrprecipG)    ! bit messy but this will internally iterate over populations
        PrprecipG=max(PrprecipG, 1e-5_wp)
      end if
      Prprecip=PrprecipG    ! glow returns rates from all populations so this becomes a straight assignments instead of accumlation
      Qeprecip=QeprecipG
      iver=iverG            ! store integrated VER computed by GLOW
    end if
  else
    !! do not compute impact ionization on a closed mesh (presumably there is no source of energetic electrons at these lats.)
    !if (mpi_cfg%myid==0 .and. debug) then
    !  print *, 'Looks like we have a closed grid, so skipping impact ionization for time step:  ',t
    !end if
  end if

  !if (mpi_cfg%myid==0) then
  !  if (debug) print *, 'Min/max root electron impact ionization production rates for time:  ',t,' :  ', &
  !    minval(Prprecip), maxval(Prprecip)
  !end if
  !if ((cfg%flagglow /= 0).and.(mpi_cfg%myid == 0)) then
  !  if (cfg%flagglow/=0 .and. debug) print *, 'Min/max 427.8 nm emission column-integrated intensity for time:  ',t,' :  ', &
  !    minval(iver(:,:,2)), maxval(iver(:,:,2))
  !end if
end subroutine impact_ionization


!> Ionization from solar radiation, *accumulates* rates, so initialize to zero if you want soley solar sources :)
subroutine solar_ionization(t,x,ymd,UTsec,f107a,f107,Prprecip,Qeprecip,ns,nn,Tn,gavg,Tninf)
  real(wp), intent(in) :: t
  class(curvmesh), intent(in) :: x
  integer, dimension(3), intent(in) :: ymd
  real(wp), intent(in) :: UTsec
  real(wp), intent(in) :: f107a,f107
  real(wp), dimension(:,:,:,:), intent(inout) :: Prprecip
  real(wp), dimension(:,:,:), intent(inout) :: Qeprecip
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: Tn
  real(wp), intent(in) :: gavg,Tninf
  real(wp), dimension(1:size(Prprecip,1),1:size(Prprecip,2),1:size(Prprecip,3),1:size(Prprecip,4)) :: Prpreciptmp
  real(wp), dimension(1:size(Qeprecip,1),1:size(Qeprecip,2),1:size(Qeprecip,3)) :: Qepreciptmp
  real(wp), dimension(1:size(Prprecip,1),1:size(Prprecip,2),1:size(Prprecip,3)) :: chi

  ! solar zenith angle
  chi=sza(ymd(1),ymd(2),ymd(3),UTsec,x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3))   ! chi size depends on glon,glat size b/c sza elemental
  !if (mpi_cfg%myid==0 .and. debug) then
  if (debug) then
    print *, 'Computing photoionization for time:  ',t,' using sza range of (root only):  ', &
      minval(chi)*180/pi, maxval(chi)*180/pi
  end if

  ! solar fluxes and resulting ionization rates
  Prpreciptmp=photoionization(x,nn,chi,f107,f107a,gavg,Tninf)
  !if (mpi_cfg%myid==0 .and. debug) then
  if (debug) then
    print *, 'Min/max root photoionization production rates for time:  ',t,' :  ', &
      minval(Prpreciptmp), maxval(Prpreciptmp)
  end if

  Prpreciptmp = max(Prpreciptmp, 1e-5_wp)
  !! enforce minimum production rate to preserve conditioning for species that rely on constant production
  !! testing should probably be done to see what the best choice is...

  Qepreciptmp = eheating(nn,Prpreciptmp,ns)
  !! thermal electron heating rate from Swartz and Nisbet, (1978)

  !> photoion ionrate and heating calculated separately, added together with ionrate and heating from Fang or GLOW
  Prprecip = Prprecip + Prpreciptmp
  Qeprecip = Qeprecip + Qepreciptmp
end subroutine solar_ionization


!> Energy source/loss solutions.  Upon entry the energy density should have the most recently updated state.  Upon exit
!    both the energy density and temperature are fully updated.
subroutine energy_source_loss(dt,Pr,Lo,Qeprecip,rhoes,Ts,ns)
  real(wp), intent(in) :: dt
  real(wp), dimension(:,:,:,:), intent(inout) :: Pr
  real(wp), dimension(:,:,:,:), intent(in) :: Lo
  real(wp), dimension(:,:,:), intent(in) :: Qeprecip
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: rhoes,Ts,ns
  real(wp), dimension(1:size(rhoes,1)-4,1:size(rhoes,2)-4,1:size(rhoes,3)-4) :: paramtrim
  integer :: isp,lsp

  lsp=size(rhoes,4)
  do isp=1,lsp
    if (isp==lsp) then
      Pr(:,:,:,lsp)=Pr(:,:,:,lsp)+Qeprecip
    end if
    paramtrim=rhoes(1:lx1,1:lx2,1:lx3,isp)
    paramtrim=ETD_uncoupled(paramtrim,Pr(:,:,:,isp),Lo(:,:,:,isp),dt)
    rhoes(1:lx1,1:lx2,1:lx3,isp)=paramtrim

    Ts(:,:,:,isp)=(gammas(isp) - 1)/kB*rhoes(:,:,:,isp)/max(ns(:,:,:,isp),mindensdiv)
    Ts(:,:,:,isp)=max(Ts(:,:,:,isp), 100._wp)
  end do
end subroutine energy_source_loss


!>  Momentum source/loss processes.  Upon entry the momentum density should be updated to most recent; upon exit
!     both momentum density and velocity will be updated.
subroutine momentum_source_loss(dt,x,Pr,Lo,ns,rhovs1,vs1)
  real(wp), intent(in) :: dt
  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:,:), intent(in) :: Pr
  real(wp), dimension(:,:,:,:), intent(in) :: Lo
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: rhovs1,vs1
  real(wp), dimension(1:size(rhovs1,1)-4,1:size(rhovs1,2)-4,1:size(rhovs1,3)-4) :: paramtrim
  real(wp), dimension(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2) :: chrgflux
  integer :: isp,lsp

!   if (maxval(abs(vs1(1:lx1,1:lx2,1:lx3,:))) > 1e4) then
!     print*, maxloc(abs(vs1(1:lx1,1:lx2,1:lx3,:)))
!     print*, 'Data corrupted before momentum source solve!'
!     print*, vs1(1:lx1,lx2,lx3,6)
!     print*, 'Data corrupted before momentum source solve!'
!     print*, vs1(1:lx1,lx2-1,lx3-1,6)
!     print*, 'Data corrupted before momentum source solve!'
!     print*, vs1(1:lx1,lx2-2,lx3-2,6)
!     error stop
!   end if

  lsp=size(rhovs1,4)
  do isp=1,lsp-1
    paramtrim=rhovs1(1:lx1,1:lx2,1:lx3,isp)
    paramtrim=ETD_uncoupled(paramtrim,Pr(:,:,:,isp),Lo(:,:,:,isp),dt)
    rhovs1(1:lx1,1:lx2,1:lx3,isp)=paramtrim
    vs1(:,:,:,isp)=rhovs1(:,:,:,isp)/(ms(isp)*max(ns(:,:,:,isp),mindensdiv))
  end do

   !print*, 'vs1 source/loss middle:  ',shape(vs1(1:lx1,1:lx2,1:lx3,:)),minval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       maxval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       minloc(vs1(1:lx1,1:lx2,1:lx3,:)),maxloc(vs1(1:lx1,1:lx2,1:lx3,:))

  ! Update velocity and momentum for electrons
  ! in keeping with the way the above situations have been handled keep the ghost cells with this calculation
  chrgflux = 0.0
  do isp=1,lsp-1
    chrgflux=chrgflux+ns(:,:,:,isp)*qs(isp)*vs1(:,:,:,isp)
  end do
  !  vs1(1:lx1,1:lx2,1:lx3,lsp)=1/max(ns(1:lx1,1:lx2,1:lx3,lsp),mindensdiv)/qs(lsp)*(J1-chrgflux)   !density floor needed???
  vs1(:,:,:,lsp)=-1/max(ns(:,:,:,lsp),mindensdiv)/qs(lsp)*chrgflux    !don't bother with FAC contribution...
  rhovs1(:,:,:,lsp)=ns(:,:,:,lsp)*ms(lsp)*vs1(:,:,:,lsp)              ! update electron momentum in case it is ever used

   !print*, 'vs1 source/loss end:  ',shape(vs1(1:lx1,1:lx2,1:lx3,:)),minval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       maxval(vs1(1:lx1,1:lx2,1:lx3,:)), &
   !                       minloc(vs1(1:lx1,1:lx2,1:lx3,:)),maxloc(vs1(1:lx1,1:lx2,1:lx3,:))
end subroutine momentum_source_loss


!> Mass source and loss processes
subroutine mass_source_loss(dt,Pr,Lo,Prprecip,ns)
  real(wp), intent(in) :: dt
  real(wp), dimension(:,:,:,:), intent(inout) :: Pr
  real(wp), dimension(:,:,:,:), intent(in) :: Lo
  real(wp), dimension(:,:,:,:), intent(in) :: Prprecip
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: paramtrim
  integer :: isp,lsp

  lsp=size(ns,4)
  Pr(:,:,:,1:6)=Pr(:,:,:,1:6)+Prprecip
  do isp=1,lsp-1
    paramtrim=ns(1:lx1,1:lx2,1:lx3,isp)
    paramtrim=ETD_uncoupled(paramtrim,Pr(:,:,:,isp),Lo(:,:,:,isp),dt)
    ns(1:lx1,1:lx2,1:lx3,isp)=paramtrim    !should there be a density floor here???  I think so...
  end do
  ns(:,:,:,lsp)=sum(ns(:,:,:,1:lsp-1),4)
end subroutine mass_source_loss


!> Deal with cells outside computation domain; i.e. apply fill values.
subroutine clean_param(x,paramflag,param)
  !------------------------------------------------------------
  !-------THIS SUBROUTINE ZEROS OUT ALL NULL CELLS AND HANDLES
  !-------POSSIBLE NULL ARTIFACTS AT BOUNDARIES
  !------------------------------------------------------------
  class(curvmesh), intent(in) :: x
  integer, intent(in) :: paramflag
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: param     !note that this is 4D and is meant to include ghost cells
  integer :: isp,ix1,ix2,ix3,iinull,ix1beg,ix1end

  select case (paramflag)
    case (1)    !density
      param(:,:,:,1:lsp-1)=max(param(:,:,:,1:lsp-1),mindens)    ! enforce a minimum density
      param(:,:,:,lsp)=sum(param(:,:,:,1:lsp-1),4)              !enforce charge neutrality based on ion densities

      do isp=1,lsp             !set null cells to some value
        if (isp==1) then
          do iinull=1,x%lnull
            ix1=x%inull(iinull,1)
            ix2=x%inull(iinull,2)
            ix3=x%inull(iinull,3)

            param(ix1,ix2,ix3,isp)=mindensnull*1e-2_wp
          end do
        else
          do iinull=1,x%lnull
            ix1=x%inull(iinull,1)
            ix2=x%inull(iinull,2)
            ix3=x%inull(iinull,3)

            param(ix1,ix2,ix3,isp)=mindensnull
          end do
        end if
      end do


      !SET DENSITY TO SOME HARMLESS VALUE in the ghost cells
      param(-1:0,:,:,:)=mindensdiv
      param(lx1+1:lx1+2,:,:,:)=mindensdiv
      param(:,-1:0,:,:)=mindensdiv
      param(:,lx2+1:lx2+2,:,:)=mindensdiv
      param(:,:,-1:0,:)=mindensdiv
      param(:,:,lx3+1:lx3+2,:)=mindensdiv

!      param(-1:0,:,:,:)=1
!      param(lx1+1:lx1+2,:,:,:)=1
!      param(:,-1:0,:,:)=1
!      param(:,lx2+1:lx2+2,:,:)=1
!      param(:,:,-1:0,:)=1
!      param(:,:,lx3+1:lx3+2,:)=1
    case (2)    !velocity
      do isp=1,lsp       !set null cells to zero mometnum
        do iinull=1,x%lnull
          ix1=x%inull(iinull,1)
          ix2=x%inull(iinull,2)
          ix3=x%inull(iinull,3)

          param(ix1,ix2,ix3,isp) = 0
        end do
      end do

      !FORCE THE BORDER CELLS TO BE SAME AS THE FIRST INTERIOR CELL (deals with some issues on dipole grids), skip for non-dipole.
      if (x%gridflag==0) then      ! closed dipole
        do isp=1,lsp
          do ix3=1,lx3
            do ix2=1,lx2
              ix1beg=1
              do while( (.not. x%nullpts(ix1beg,ix2,ix3)) .and. ix1beg<lx1)     !find the first non-null index for this field line, need to be careful if no null points exist...
                ix1beg=ix1beg+1
              end do

              ix1end=ix1beg
              do while(x%nullpts(ix1end,ix2,ix3) .and. ix1end<lx1)     !find the last non-null index for this field line
                ix1end=ix1end+1
              end do
              !if (ix1end /= ix1beg .and. ix1end /= lx1) ix1end=ix1end-1      ! I think this has been left out for a long time!?

              if (ix1beg /= lx1) then    !only do this if we actually have null grid points
                param(ix1beg,ix2,ix3,isp)=param(ix1beg+1,ix2,ix3,isp)
              end if
              if (ix1end /= lx1) then
                param(ix1end,ix2,ix3,isp)=param(ix1end-1,ix2,ix3,isp)
              end if
            end do
          end do
        end do
      elseif (x%gridflag==1) then     ! open dipole grid, inverted
        do isp=1,lsp
          do ix3=1,lx3
            do ix2=1,lx2
              ix1end=1
              do while((.not. x%nullpts(ix1end,ix2,ix3)) .and. ix1end<lx1)     !find the first non-null index for this field line
                ix1end=ix1end+1
              end do

              if (ix1end /= lx1) then
                param(ix1end,ix2,ix3,isp)=param(ix1end-1,ix2,ix3,isp)
              end if
            end do
          end do
        end do
      end if

  !MZ - for reasons I don't understand, this causes ctest to fail...  Generates segfaults everywhere in the CI (these are due to failing the comparisons)...  Okay so the deal here is that the ghost cell velocity values are used to compute artificial viscosity in fluid_adv, so one cannot clear them out without ruining the solution.  AFAIK no other params have this issue...
      !ZERO OUT THE GHOST CELL VELOCITIES
  !    param(-1:0,:,:,:)= 0
  !    param(lx1+1:lx1+2,:,:,:)= 0
  !    param(:,-1:0,:,:)= 0
  !    param(:,lx2+1:lx2+2,:,:)= 0
  !    param(:,:,-1:0,:)= 0
  !    param(:,:,lx3+1:lx3+2,:)= 0
    case (3)    !temperature
      param=max(param,100._wp)     !temperature floor

      do isp=1,lsp       !set null cells to some value
        do iinull=1,x%lnull
          ix1=x%inull(iinull,1)
          ix2=x%inull(iinull,2)
          ix3=x%inull(iinull,3)

          param(ix1,ix2,ix3,isp) = 100
        end do
      end do

      !> SET TEMPS TO SOME NOMINAL VALUE in the ghost cells
      param(-1:0,:,:,:) = 100
      param(lx1+1:lx1+2,:,:,:) = 100
      param(:,-1:0,:,:) = 100
      param(:,lx2+1:lx2+2,:,:) = 100
      param(:,:,-1:0,:) = 100
      param(:,:,lx3+1:lx3+2,:) = 100
    case default
      !! throw an error as the code is likely not going to behave in a predictable way in this situation...
      error stop '!non-standard parameter selected in clean_params, unreliable/incorrect results possible...'
  end select
end subroutine clean_param


!> Deal with cells outside computation domain; do not touch ghost cells in any way - exercise caution in the way null cells are
!    treated as compared to a "normal" clean.  Screen cells "near" null cells for excessively large parameter values that may
!    result from interpolation artifacts.
subroutine clean_param_after_regrid(x,paramflag,param,Tn)
  !------------------------------------------------------------
  !-------THIS SUBROUTINE ZEROS OUT ALL NULL CELLS AND HANDLES
  !-------POSSIBLE NULL ARTIFACTS AT BOUNDARIES
  !------------------------------------------------------------
  class(curvmesh), intent(in) :: x
  integer, intent(in) :: paramflag
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: param     !note that this is 4D and is meant to include ghost cells
  real(wp), dimension(:,:,:), intent(in) :: Tn
  integer :: isp,ix1,ix2,ix3,iinull,ix1beg,ix1end,ix2beg
  integer :: ibuf
  integer, parameter :: lbuf=3     
  !^ controls how far from null cells we want to forcibly replace data after a refine.  What is effectively being done here
  !    is that we are copying data from some specified location into the buffer region.  This is likely acceptable because
  !    the regions where we do this are very close to equilibrium so they will "snap" back to that state rather than generating
  !    some weird transient from the forcibly replacement.  

  select case (paramflag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Density
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (1)
      param(:,:,:,1:lsp-1)=max(param(:,:,:,1:lsp-1),mindens)    ! enforce a minimum density
      param(:,:,:,lsp)=sum(param(:,:,:,1:lsp-1),4)              !enforce charge neutrality based on ion densities

      do isp=1,lsp             !set null cells to some value
        if (isp==1) then
          do iinull=1,x%lnull
            ix1=x%inull(iinull,1)
            ix2=x%inull(iinull,2)
            ix3=x%inull(iinull,3)

            param(ix1,ix2,ix3,isp)=mindensnull*1e-2_wp
          end do
        else
          do iinull=1,x%lnull
            ix1=x%inull(iinull,1)
            ix2=x%inull(iinull,2)
            ix3=x%inull(iinull,3)

            param(ix1,ix2,ix3,isp)=mindensnull
          end do
        end if
      end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Drift velocity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (2)
      do isp=1,lsp       !set null cells to zero mometnum
        do iinull=1,x%lnull
          ix1=x%inull(iinull,1)
          ix2=x%inull(iinull,2)
          ix3=x%inull(iinull,3)

          param(ix1,ix2,ix3,isp) = 0
        end do
      end do

      !FORCE THE BORDER CELLS TO BE SAME AS THE FIRST INTERIOR CELL (deals with some issues on dipole grids), skip for non-dipole.
      if (x%gridflag==0) then      ! closed dipole
        do isp=1,lsp
          do ix3=1,lx3
            do ix2=1,lx2
              ix1beg=1
              do while( (.not. x%nullpts(ix1beg,ix2,ix3)) .and. ix1beg<lx1)     !find the first non-null index for this field line, need to be careful if no null points exist...
                ix1beg=ix1beg+1
              end do

              ix1end=ix1beg
              do while(x%nullpts(ix1end,ix2,ix3) .and. ix1end<lx1)     !find the last non-null index for this field line
                ix1end=ix1end+1
              end do

              if (ix1beg /= lx1) then    !only do this if we actually have null grid points
                do ibuf=1,lbuf
                  !param(ix1beg+ibuf-1,ix2,ix3,isp)=param(ix1beg+lbuf,ix2,ix3,isp)
                  param(ix1beg+ibuf-1,ix2,ix3,isp)=0._wp
                end do
              end if
              if (ix1end /= lx1) then
                do ibuf=1,lbuf
                  !param(ix1end+ibuf-1,ix2,ix3,isp)=param(ix1end-lbuf,ix2,ix3,isp)
                  param(ix1end+ibuf-1,ix2,ix3,isp)=0._wp
                end do
              end if
            end do
          end do
        end do
      elseif (x%gridflag==1) then     ! open dipole grid, inverted
        do isp=1,lsp
          do ix3=1,lx3
            do ix2=1,lx2
              ix1end=1
              do while((.not. x%nullpts(ix1end,ix2,ix3)) .and. ix1end<lx1)     !find the first non-null index for this field line
                ix1end=ix1end+1
              end do

              if (ix1end /= lx1) then
                param(ix1end,ix2,ix3,isp)=param(ix1end-1,ix2,ix3,isp)
              end if
            end do
          end do
        end do
      end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Temperature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (3)
      param=max(param,100._wp)     !temperature floor

      do isp=1,lsp       !set null cells to some value
        do iinull=1,x%lnull
          ix1=x%inull(iinull,1)
          ix2=x%inull(iinull,2)
          ix3=x%inull(iinull,3)

          param(ix1,ix2,ix3,isp) = 100
        end do
      end do


      !FORCE THE BORDER CELLS TO BE SAME AS THE FIRST INTERIOR CELL (deals with some issues on dipole grids), skip for non-dipole.
      if (x%gridflag==0) then      ! closed dipole
        do isp=1,lsp
          do ix3=1,lx3
            do ix2=1,lx2
              ix1beg=1
              do while( (.not. x%nullpts(ix1beg,ix2,ix3)) .and. ix1beg<lx1)     !find the first non-null index for this field line, need to be careful if no null points exist...
                ix1beg=ix1beg+1
              end do

              ix1end=ix1beg
              do while(x%nullpts(ix1end,ix2,ix3) .and. ix1end<lx1)     !find the last non-null index for this field line
                ix1end=ix1end+1
              end do

              if (ix1beg /= lx1) then    !only do this if we actually have null grid points
                do ibuf=1,lbuf
                  !param(ix1beg+ibuf-1,ix2,ix3,isp)=param(ix1beg+lbuf,ix2,ix3,isp)
                  param(ix1beg+ibuf-1,ix2,ix3,isp)=Tn(ix1beg+ibuf-1,ix2,ix3)
                end do
              end if
              if (ix1end /= lx1) then
                do ibuf=1,lbuf
                  !param(ix1end+ibuf-1,ix2,ix3,isp)=param(ix1end-lbuf,ix2,ix3,isp)
                  param(ix1end+ibuf-1,ix2,ix3,isp)=Tn(ix1end+ibuf-1,ix2,ix3)
                end do
              end if
            end do
          end do
        end do
      elseif (x%gridflag==1) then     ! open dipole grid, inverted
        do isp=1,lsp
          do ix3=1,lx3
            do ix2=1,lx2
              ix1end=1
              do while((.not. x%nullpts(ix1end,ix2,ix3)) .and. ix1end<lx1)     !find the first non-null index for this field line
                ix1end=ix1end+1
              end do

              if (ix1end /= lx1) then
                param(ix1end,ix2,ix3,isp)=param(ix1end-1,ix2,ix3,isp)
              end if
            end do
          end do
        end do
      end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Error
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case default
      !! throw an error as the code is likely not going to behave in a predictable way in this situation...
      error stop '!non-standard parameter selected in clean_params, unreliable/incorrect results possible...'
  end select
end subroutine clean_param_after_regrid

end module multifluid
