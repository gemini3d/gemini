module collisions

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp, lsp, ln, ms, mn, kb, pi, elchrg, qs, debug
use gemini3d_config, only: gemini_cfg
use meshobj, only : curvmesh


implicit none (type, external)
private
public :: thermal_conduct, conductivities, capacitance, maxwell_colln, coulomb_colln, NLConductivity


real(wp), parameter :: Csn(lsp,ln) = reshape( [ real(wp) :: &
-1, 6.82e-10, 6.64e-10, -1, &
2.44e-10, 4.34e-10, 4.27e-10, 0.69e-10, &
2.58e-10, -1, 4.49e-10, 0.74e-10, &
2.31e-10, 4.13e-10, -1, 0.65e-10, &
4.42e-10, 7.47e-10, 7.25e-10, 1.45e-10, &
-1, 33.6e-10, 32.0e-10, -1, &
-1, -1, -1, -1], shape(Csn), order=[2,1])

real(wp), parameter :: C2sn1(lsp,ln) = reshape( &
[real(wp) :: &
3.67e-11, 0, 0, 4.63e-12, &
0, 0, 0, 0, &
0, 5.14e-11, 0, 0, &
0, 0, 2.59e-11, 0, &
0, 0, 0, 0, &
6.61e-11, 0, 0, 2.65e-10, &
-1, -1, -1, -1], shape(C2sn1), order=[2,1])

real(wp), parameter :: C2sn2(lsp,ln) = reshape( &
[real(wp) :: &
0.064, 0, 0, -1, &
0, 0, 0, 0, &
0, 0.069, 0, 0, &
0, 0, 0.073, 0, &
0, 0, 0, 0, &
0.047, 0, 0, 0.083, &
-1, -1, -1, -1], shape(C2sn2), order=[2,1])

real(wp), parameter :: Csj(lsp,lsp) = reshape( &
[real(wp) :: &
0.22, 0.26, 0.25, 0.26, 0.22, 0.077, 1.87e-3, &
0.14, 0.16, 0.16, 0.17, 0.13, 0.042, 9.97e-4, &
0.15, 0.17, 0.17, 0.18, 0.14, 0.045, 1.07e-3, &
0.13, 0.16, 0.15, 0.16, 0.12, 0.039, 9.347e-4, &
0.25, 0.28, 0.28, 0.28, 0.24, 0.088, 2.136e-3, &
1.23, 1.25, 1.25, 1.25, 1.23, 0.90,  29.7e-3, &
54.5, 54.5, 54.5, 54.5, 54.5, 54.5,  38.537], shape(Csj), order=[2,1])

real(wp), parameter :: thermal_coeff(lsp-1) = [real(wp) :: &
0.1019e-12, 0.0747e-12, 0.0754e-12, 0.0701e-12, &
0.1068e-12, 0.3986e-12]

contains


subroutine maxwell_colln(isp,isp2,nn,Tn,Ts,nusn)
  !! COMPUTE MAXWELL COLLISIONS OF ISP WITH ISP2.  ION
  !! TEMPERATURE/DENSITY ARRAYS EXPECTED TO INCLUDE GHOST CELLS
  !! Note that it is done on a per species basis
  integer, intent(in) :: isp,isp2
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: Tn
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: Ts
    real(wp), dimension(1:size(Tn,1),1:size(Tn,2),1:size(Tn,3)), intent(inout) :: nusn
  !! intent(out)
  integer :: lx1,lx2,lx3
  real(wp),dimension(1:size(Tn,1),1:size(Tn,2),1:size(Tn,3)) :: Teff
  
  lx1=size(Ts,1)-4
  lx2=size(Ts,2)-4
  lx3=size(Ts,3)-4
  
  
  if (isp<lsp) then
  !! ion-neutral
    if (Csn(isp,isp2) < 0) then
    !! resonant
      if (isp==1 .and. isp2==4) then
        Teff=Tn+Ts(1:lx1,1:lx2,1:lx3,isp) / 16
        nusn=C2sn1(isp,isp2)*Teff**0.5*nn(:,:,:,isp2)*1e-6_wp
      else
        Teff=0.5*(Tn+Ts(1:lx1,1:lx2,1:lx3,isp))
        nusn=C2sn1(isp,isp2)*(1 - C2sn2(isp,isp2)*log10(Teff))**2 * &
             (Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
      end if
    else
    !! nonresonant
      nusn=Csn(isp,isp2)*nn(:,:,:,isp2)*1e-6_wp
    end if
  else
  ! !! electron-neutral
  !   Teff=Ts(1:lx1,1:lx2,1:lx3,isp)
  
  !   select case (isp2)
  !     case (1)
  !       nusn=8.9e-11_wp*(1.0+5.7e-4_wp*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
  !     case (2)
  !       nusn=2.33e-11_wp*(1.0-1.21e-4_wp*Teff)*(Teff)*nn(:,:,:,isp2)*1e-6_wp
  !     case (3)
  !       nusn=1.82e-10_wp*(1.0+3.6e-2_wp*(Teff**0.5))*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
  !     case (4)
  !       nusn=4.5e-9_wp*(1.0-1.35e-4_wp*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
  !     case default
  !       write(stderr,*) 'ERROR: isp2 value is unknown: ',isp2
  !       error stop
  !   end select
  ! end if

!!Matt fix, caps temperature when calculating cross sections (?)
!! electron-neutral
    Teff=Ts(1:lx1,1:lx2,1:lx3,isp)
    select case (isp2)
      case (1)
        nusn=8.9e-11_wp*(1.0+5.7e-4_wp*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
      case (2)
        where (Teff>5000._wp)    ! avoids negative collision frequency!
          Teff=5000._wp
        end where

        nusn=2.33e-11_wp*(1.0-1.21e-4_wp*Teff)*(Teff)*nn(:,:,:,isp2)*1e-6_wp
      case (3)
        nusn=1.82e-10_wp*(1.0+3.6e-2_wp*(Teff**0.5))*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
      case (4)
        where (Teff>5000._wp)
          Teff=5000._wp
        end where

        nusn=4.5e-9_wp*(1.0-1.35e-4_wp*Teff)*(Teff**0.5)*nn(:,:,:,isp2)*1e-6_wp
      case default
        write(stderr,*) 'ERROR: isp2 value is unknown: ',isp2
        error stop
    end select
  end if
end subroutine maxwell_colln


pure subroutine coulomb_colln(isp,isp2,ns,Ts,vs1,nusj,Phisj,Psisj)
  !! COMPUTE COULOMB COLLISIONS OF ISP WITH ISP2.
  !! TEMPERATURE/DENSITY ARRAYS EXPECTED TO INCLUDE GHOST CELLS
  !! NOTE THAT OTHER PIECES OF THE CODE REQUIRE SELF COLLISIONS
  !! TO BE ZERO TO YIELD CORRECT OUTPUT (SOURCES.MOD)
  !! Note that it is done on a per species basis
  
  integer, intent(in) :: isp,isp2
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
  
  real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), intent(inout) :: nusj,Phisj,Psisj
  !! intent(out)
  integer :: lx1,lx2,lx3
  real(wp) :: mred
  real(wp),dimension(1:size(Ts,1)-4, 1:size(Ts,2)-4, 1:size(Ts,3)-4) &
            :: Teff,Wsj
  
  
  lx1=size(Ts,1)-4
  lx2=size(Ts,2)-4
  lx3=size(Ts,3)-4
  
  if (isp==isp2) then
  !! zero out all self collision terms (would need to be changed if non-Maxwellian distribution used).
    nusj = 0
    Phisj = 0
    Psisj = 0
  else
    Teff=(ms(isp2)*Ts(1:lx1,1:lx2,1:lx3,isp)+ms(isp)* &
           Ts(1:lx1,1:lx2,1:lx3,isp2))/(ms(isp2)+ms(isp))
    nusj=Csj(isp,isp2)*ns(1:lx1,1:lx2,1:lx3,isp2)*1e-6_wp/Teff**1.5_wp


    mred=ms(isp)*ms(isp2)/(ms(isp)+ms(isp2))
    Wsj=abs(vs1(1:lx1,1:lx2,1:lx3,isp)-vs1(1:lx1,1:lx2,1:lx3,isp2))/ &
          sqrt(2*kB*Teff/mred)
    Psisj=exp(-Wsj**2)
    where (Wsj<0.1_wp)
      Phisj=1
    elsewhere
      Phisj=3.0_wp/4*sqrt(pi)*erf(Wsj)/Wsj**3 - 3.0_wp/2/Wsj**2*Psisj
    end where
  end if
end subroutine coulomb_colln


! pure subroutine thermal_conduct(isp,Ts,ns,nn,J1,lambda,beta)
!   !! COMPUTE THERMAL CONDUCTIVITY.
!   !! TEMPERATURE ARRAY IS EXPECTED TO INCLUDE GHOST CELLS
!   !! Note that it is done on a per species basis
  
!   integer, intent(in) :: isp
!   real(wp), dimension(-1:,-1:,-1:), intent(in) :: Ts,ns
  
!   real(wp), dimension(:,:,:,:), intent(in) :: nn
!   real(wp), dimension(-1:,-1:,-1:), intent(in) :: J1
  
!   real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), intent(inout) :: lambda,beta
!   !! intent(out)
  
!   integer :: lx1,lx2,lx3
  
  
!   lx1=size(Ts,1)-4
!   lx2=size(Ts,2)-4
!   lx3=size(Ts,3)-4
  
!   if (isp<lsp) then
!   !! ion species
!   !  lambda=25.0_wp/8 * kB**2*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2)/ms(isp)/(Csj(isp,isp)*1e-6_wp)
!     !! avoids precision issues by precomputing the transport coefficients (see parameter blocks above)
!     lambda=thermal_coeff(isp)*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2)
!     beta=0.0
!   else                  !electrons
!     lambda=elchrg * 100 * 7.7e5_wp*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2)/ &
!         (1 + 3.22e4_wp*Ts(1:lx1,1:lx2,1:lx3)**2/ns(1:lx1,1:lx2,1:lx3)* &
!           (nn(:,:,:,1)*1.1e-16_wp*(1+5.7e-4_wp*Ts(1:lx1,1:lx2,1:lx3)) + &
!           nn(:,:,:,2)*2.82e-17_wp*sqrt(Ts(1:lx1,1:lx2,1:lx3))* &
!           (1-1.21e-4_wp*Ts(1:lx1,1:lx2,1:lx3))+nn(:,:,:,3)* &
!           2.2e-16_wp*(1+3.6e-2_wp*sqrt(Ts(1:lx1,1:lx2,1:lx3))) ))
!     beta=5.0_wp/2 * kB/elchrg * J1(1:lx1,1:lx2,1:lx3)
!   end if
! end subroutine thermal_conduct

! Matt version that caps thermal conduct
pure subroutine thermal_conduct(isp,Ts,ns,nn,J1,lambda,beta)
  !! COMPUTE THERMAL CONDUCTIVITY.
  !! TEMPERATURE ARRAY IS EXPECTED TO INCLUDE GHOST CELLS
  !! Note that it is done on a per species basis
  
  integer, intent(in) :: isp
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: Ts,ns
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: J1
  real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4), intent(inout) :: lambda,beta
  !! intent(out)
  real(wp), dimension(1:size(Ts,1)-4,1:size(Ts,2)-4,1:size(Ts,3)-4) :: Tstmp
  integer :: lx1,lx2,lx3
  
  
  lx1=size(Ts,1)-4
  lx2=size(Ts,2)-4
  lx3=size(Ts,3)-4
  
  if (isp<lsp) then
  !! ion species
  !  lambda=25.0_wp/8 * kB**2*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2)/ms(isp)/(Csj(isp,isp)*1e-6_wp)
    !! avoids precision issues by precomputing the transport coefficients (see parameter blocks above)
    lambda=thermal_coeff(isp)*Ts(1:lx1,1:lx2,1:lx3)**(5.0_wp/2)
    beta=0.0
  else                  !electrons
    Tstmp=Ts(1:lx1,1:lx2,1:lx3)
    where (Ts(1:lx1,1:lx2,1:lx3) > 6000.0)
      Tstmp=6000.0
    end where
    lambda=elchrg * 100 * 7.7e5_wp*Tstmp**(5.0_wp/2)/ &
        (1 + 3.22e4_wp*Tstmp**2/ns(1:lx1,1:lx2,1:lx3)* &
          (nn(:,:,:,1)*1.1e-16_wp*(1+5.7e-4_wp*Tstmp) + &
          nn(:,:,:,2)*2.82e-17_wp*sqrt(Tstmp)* &
          (1-1.21e-4_wp*Tstmp)+nn(:,:,:,3)* &
          2.2e-16_wp*(1+3.6e-2_wp*sqrt(Tstmp) )))
    beta=5.0_wp/2 * kB/elchrg * J1(1:lx1,1:lx2,1:lx3)
  end if
end subroutine thermal_conduct


subroutine conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
  !! COMPUTE THE CONDUCTIVITIES OF THE IONOSPHERE.  STATE
  !! VARS. INCLUDE GHOST CELLS
  
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: Tn
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(inout) :: sig0,sigP,sigH
  !! intent(out)
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,lsp), intent(inout) :: muP,muH
  !! intent(out)
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,lsp), intent(inout) :: nusn
  !! intent(out)
  !! defined for each ion species, summed over neutral species
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(inout) :: sigPgrav,sigHgrav
  !! intent(out)
  
  integer :: isp,isp2,lx1,lx2,lx3
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: OMs
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: nuej,Phisj,Psisj,nutmp,mupar,mubase,rho
  
  lx1=size(Ts,1)-4
  lx2=size(Ts,2)-4
  lx3=size(Ts,3)-4
  
  !> Refactor this code so that it outputs nusn instead of the two "neutral mobilities", this also facilitates pressure terms...
  
  !MOBILITIES
  do isp=1,lsp
  !      OMs=qs(isp)*abs(B1)/ms(isp)
    !! cyclotron, abs() is sketch, needs to be checked.
    !! Basically a negative sign here is fine, while abs messes up direction of Hall current
    OMs=qs(isp)*B1(1:lx1,1:lx2,1:lx3)/ms(isp)
    !! cyclotron, a negative sign from B1 here is fine for cartesian, but for dipole this should be the magnitude
    !! since the magnetic field is *assumed* to be along the x1-direction
  
    nusn(:,:,:,isp) = 0
    do isp2=1,ln
      call maxwell_colln(isp,isp2,nn,Tn,Ts,nutmp)
      nusn(:,:,:,isp)=nusn(:,:,:,isp)+nutmp
    end do
  
    if (isp<lsp) then
      mubase=qs(isp)/ms(isp)/nusn(:,:,:,isp)      !parallel mobility
    else
      nuej = 0
      do isp2=1,lsp
        call coulomb_colln(isp,isp2,ns,Ts,vs1,nutmp,Phisj,Psisj)
        nuej=nuej+nutmp
      end do
  
      mupar=qs(lsp)/ms(lsp)/(nusn(:,:,:,isp)+nuej)
      mubase=qs(lsp)/ms(lsp)/nusn(:,:,:,isp)
    end if
  
    !! modified mobilities for neutral wind calculations.
    !! these are deprecated since we output collision freq.
  !  muPvn(:,:,:,isp)=nu**2/(nu**2+OMs**2)
  !  muHvn(:,:,:,isp)= -nu*OMs/(nu**2+OMs**2)
  
    !electrical mobilities
    muP(:,:,:,isp)=mubase*nusn(:,:,:,isp)**2/(nusn(:,:,:,isp)**2+OMs**2)                !Pederson
  !  if (isp==lsp) then
  !    print*, 'Flipping electron mobility signs:  ',lsp
  !    print*, minval(muP(:,:,:,lsp)),maxval(muP(:,:,:,lsp))
  !    muP(:,:,:,lsp)=-1.0*muP(:,:,:,lsp)
  !    print*, minval(muP(:,:,:,lsp)),maxval(muP(:,:,:,lsp))
  !  end if
    muH(:,:,:,isp) = mubase*-1*nusn(:,:,:,isp)*OMs/(nusn(:,:,:,isp)**2+OMs**2)       !Hall
  
    !gravity mobilities???
  end do
  
  
  !> CONDUCTIVITIES
  sig0=ns(1:lx1,1:lx2,1:lx3,lsp)*qs(lsp)*mupar
  !! parallel includes only electrons...
  
  sigP = 0
  sigH = 0
  do isp=1,lsp
    rho=ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)
    !! rho is charge density here
    sigP=sigP+rho*muP(:,:,:,isp)
    sigH=sigH+rho*muH(:,:,:,isp)
  end do
  
  !    sigH=max(sigH,0.0_wp)
  !! to deal with precision issues.
  !! This actually causes errors in Cartesian northern hemisphere grids...
  
  
  !Gravitational "conductivities"
  sigPgrav = 0
  sigHgrav = 0
  do isp=1,lsp
    rho=ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)
    !! here, rho is used as mass density
    sigPgrav=sigPgrav+rho*muP(:,:,:,isp)
    sigHgrav=sigHgrav+rho*muH(:,:,:,isp)
  end do
end subroutine conductivities


subroutine NLConductivity(nn,Tn,ns,Ts,E2,E3,x,sigP,sigH,sigNCP,sigNCH)
  !! Inputs Needed
  real(wp), dimension(:,:,:,:), intent(in) :: nn !Neutral density
  real(wp), dimension(:,:,:), intent(in) :: Tn !neutral temperature
  real(wp), dimension(:,:,:), intent(in) :: sigP,sigH !linear conductivities
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts !Plasma density and temperature
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: E2,E3 !Electric Field
  class(curvmesh), intent(in) :: x !Grid, doing this because BMAG is stored here, added at the top of the file too
  
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: sigNCP,sigNCH  
  
  !!Internal Arrays
  integer :: isp,isp2,lx1,lx2,lx3
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp) :: nsuAvg 
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp-1) :: niW 
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,ln) :: nuW 
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,2) :: nuAvg, msAvg, TsAvg  
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: Bmagnitude, nu, nsAvg, omegae, omegai, ki, ke, phi
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: Eth0, Ethreshold, Emagnitude, commonfactor
  integer, dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: FBIbinary
  
  
  !!Start
  lx1=x%lx1
  lx2=x%lx2
  lx3=x%lx3
  Bmagnitude=x%Bmag(1:lx1,1:lx2,1:lx3)
  Emagnitude=sqrt(E2(1:lx1,1:lx2,1:lx3)**2+E3(1:lx1,1:lx2,1:lx3)**2) !!Already evaluated with no ghost cells
  
  !!Initialize arrays as 0s
  nuAvg=0.0_wp
  nsuAvg=0.0_wp
  msAvg=0.0_wp
  nsAvg=0.0_wp
  TsAvg=0.0_wp
  sigNCH=0.0_wp
  sigNCP=0.0_wp
  FBIbinary=1
  
  
  !MassDensity Weight of Neutrals
  do isp2=1,ln
      nuW(:,:,:,isp2)=nn(:,:,:,isp2)*mn(isp2) !Weight of the neutrals
  end do
  
  !!MassDensity Weight of all ions
  do isp=1,lsp-1
    niW(:,:,:,isp)=ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)
  end do
  
  !! Average Collisuons frequencies: first averaging over neutrals
  do isp=1,lsp
    do isp2=1,ln
      call maxwell_colln(isp,isp2,nn,Tn,Ts,nu)
      nsuAvg(:,:,:,isp)=nsuAvg(:,:,:,isp)+nu*nuW(:,:,:,isp2) !Store the collision frequencies weighted by massdensity
    end do
    nsuAvg(:,:,:,isp)=nsuAvg(:,:,:,isp)/sum(nuW, dim=4) !! Average over all neutrals weighted by MassDensity
  end do
  
  !Final ion neutral collision frequency using only NO+ and O2+
  !Store summation of collision frequencies weighted by MassDensity
  nuAvg(:,:,:,1)=(nsuAvg(:,:,:,2)*niW(:,:,:,2)+nsuAvg(:,:,:,4)*niW(:,:,:,4))/(niW(:,:,:,2)+niW(:,:,:,4))
  
  !!Electrons do not need averaging
  nuAvg(:,:,:,2)=nsuAvg(:,:,:,lsp)
  
  !! Average mass of ions, also weighted by MassDensity
  msAvg(:,:,:,1)=(ms(2)*niW(:,:,:,2)+ms(4)*niW(:,:,:,4))/(niW(:,:,:,2)+niW(:,:,:,4))
  
  msAvg(:,:,:,2)=ms(lsp) !! Electron mass
  
  !! Average density
  !! Average just O2+ and NO+
  nsAvg=(ns(1:lx1,1:lx2,1:lx3,2)*niW(:,:,:,2)+ns(1:lx1,1:lx2,1:lx3,4)*niW(:,:,:,4))/(niW(:,:,:,2)+niW(:,:,:,4))
  
  !! ki value
  omegai=elchrg*Bmagnitude/msAvg(:,:,:,1) !! Would this work?, it will, I defined Bmagnitude above
  ki=abs(omegai/nuAvg(:,:,:,1)) !! Could do ABS to be sure of the sign
  !! ke value
  omegae=elchrg*Bmagnitude/msAvg(:,:,:,2) 
  ke=abs(omegae/nuAvg(:,:,:,2)) !!Not sure anymore about the ABS, have to ask Meers
  
  !!Phi value 1/(ki*ke)
  phi=1.0_wp/(ke*ki)
  
  !!Average ion temperature
  TsAvg(:,:,:,1)=(Ts(1:lx1,1:lx2,1:lx3,2)*niW(:,:,:,2)+Ts(1:lx1,1:lx2,1:lx3,4)*niW(:,:,:,4))/(niW(:,:,:,2)+niW(:,:,:,4))
  TsAvg(:,:,:,2)=Ts(1:lx1,1:lx2,1:lx3,lsp)
  
  !!Ethreshold
  !Ethresholdnum=(1+phi)*Bmagnitude*SQRT(kB*(1+ki**2)*(TsAvg(:,:,:,1)+TsAvg(:,:,:,2)))
  !Ethresholdden=SQRT((1-ki**2)*msAvg(:,:,:,1)) 
  
  !doi:10.1029/2011JA016649
  Eth0=20.0_wp*SQRT((TsAvg(:,:,:,1)+TsAvg(:,:,:,2))/600.0_wp)*(Bmagnitude/5.0e-5_wp) !B is written as 5e4nT, to T
  Ethreshold=(1.0_wp+phi)*SQRT((1.0_wp+ki**2)/(1.0_wp-ki**2))*Eth0*1.0e-3_wp !the 1e-3 is needed since this eq gives mV/m, not V/m
  
  !Create matrix of 1 and 0s where FBI is possible, FBIbinary starts with all 1's meaning FBI everywhere
  where (Emagnitude<=Ethreshold) !Anything without a sufficiente E field gets back to normal.
    FBIbinary=0
  end where
  
  where (ki>1.0_wp) !Anything where ions are magnetized also goes back to normal
    FBIbinary=0
  end where
  
  !Calculate conductivity term only where FBI is possible
  where (FBIbinary==1)
    commonfactor=(Emagnitude/Ethreshold-1)*(1-Ethreshold/Emagnitude)
    sigNCP=(1-ki**2)*commonfactor*sigP/(1+ki**2)
    sigNCH=-(2*ki)*(1+phi)*commonfactor*sigH/(1+ki**2)
  end where
end subroutine NLConductivity


subroutine capacitance(ns,B1,cfg,incap)
  !! COMPUTE THE INERTIAL CAPACITANCE OF THE IONOSPHERE.
  !! DENSITY/MAG FIELD STATE VARIABLE INCLUDES GHOST CELLS.
  
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
  type(gemini_cfg), intent(in) :: cfg
  real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), intent(inout) :: incap
  !! intent(out)
  integer :: lx1,lx2,lx3,isp
  
  lx1=size(ns,1)-4
  lx2=size(ns,2)-4
  lx3=size(ns,3)-4
  
  incap = 0
  do isp=1,lsp
    incap=incap+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)
  end do
  
  incap=incap/B1(1:lx1,1:lx2,1:lx3)**2
  
  if (cfg%flagcap==2) then
    if (debug) print *, '!!! Augmenting capacitance with a magnetospheric contribution...'
    incap=incap + cfg%magcap / 980e3_wp
    !! augment the value to account for a magnetosheric contribution, based on user input.  Probably should
    !! be in the 5-35 F range...  Note that the assumes that the grid extends from ~90-1000 km in altitude approx.
  end if
end subroutine capacitance


end module collisions
