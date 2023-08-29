module sources

use calculus, only : grad3d1
use collisions, only:  maxwell_colln, coulomb_colln
use phys_consts, only: wp, lsp, amu, kb, qs, ln, ms, gammas, elchrg, mn
use meshobj, only : curvmesh
use gemini3d_config, only: gemini_cfg


implicit none (type, external)
private
public :: srcsenergy, srcsmomentum, srcscontinuity

interface srcsMomentum
  module procedure srcsMomentum_curv
end interface srcsMomentum

contains


pure subroutine srcsContinuity(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts, Pr, Lo)

!------------------------------------------------------------
!-------POPULATE SOURCE/LOSS ARRAYS FOR CONTINUITY EQUATION.  ION
!-------PARAMETER ARGUMENTS (AND GRID STUFF) SHOULD INCLUDE GHOST CELLS
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(inout) :: Pr,Lo
!! intent(out)

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: betanow,kreac,Teff,Te,dv2
integer :: lx1,lx2,lx3

lx1=size(ns,1)-4
lx2=size(ns,2)-4
lx3=size(ns,3)-4

Pr=0
Lo=0
Te=Ts(1:lx1,1:lx2,1:lx3,lsp)    !< Used in calculation of Lo
dv2=(vs1(1:lx1,1:lx2,1:lx3,1)-vn1)**2+(vs2(1:lx1,1:lx2,1:lx3,1)-vn2)**2+ &
     (vs3(1:lx1,1:lx2,1:lx3,1)-vn3)**2    !gets used several times in this subprogram



!!!!!!!!!!!!!!!!!!!!!!!!!!! O+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O+ + N2 --> NO+ + N
Teff=28/(16+28._wp)*(16*amu/3/kB*(dv2) &
     + Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn


where (Teff<=3725)
  kreac=1.71676e-12_wp &
    -7.19934e-13_wp*(Teff/300) &
    +1.33276e-13_wp*(Teff/300)**2 &
    -9.28213e-15_wp*(Teff/300)**3 &
    +6.39557e-16_wp*(Teff/300)**4
end where
where (Teff>3725 .and. Teff<=30000)
  kreac=-1.52489e-11_wp &
    +7.67112e-13_wp*(Teff/300) &
    +1.19064e-13_wp*(Teff/300)**2 &
    -1.30858e-15_wp*(Teff/300)**3 &
    +4.67756e-18_wp*(Teff/300)**4
end where
!If it is above 70.000K we are doing something wrong in my opinion (Joaquin)
where (Teff>30000 .and. Teff<=70000)
  kreac=-3.2999846-9_wp &
    +3.7832649e-13_wp*(Teff) &
    -1.5807103e-17_wp*(Teff)**2 &
    +3.5017809e-22_wp*(Teff)**3 &
    -4.3053426e-27_wp*(Teff)**4 &
    +2.7760068e-32_wp*(Teff)**5 &
    -7.3160029e-38_wp*(Teff)**6 
end where

! where (Teff>30000)
!   kreac=-1.52489e-11_wp &
!     +7.67112e-13_wp*(100) &
!     +1.19064e-13_wp*(100)**2 &
!     -1.30858e-15_wp*(100)**3 &
!     +4.67756e-18_wp*(100)**4
! end where



betanow=kreac*nn(:,:,:,2)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + O2 --> O2+ + O
Teff=32/(16+32._wp)*(16*amu/3/kB*(dv2) &
     +Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn

where (Teff<=4800)
  kreac=2.78932e-11_wp &
    -6.92612e-12_wp*(Teff/300) &
    +8.67684e-13_wp*(Teff/300)**2 &
    -3.47251e-14_wp*(Teff/300)**3 &
    +5.07097e-16_wp*(Teff/300)**4
end where
where (Teff>4800 .and. Teff<=30000)
  kreac=-1.74046e-11_wp &
    +3.02328e-12_wp*Teff/300 &
    -2.39214e-15_wp*(Teff/300)**2 &
    -4.02394e-17_wp*(Teff/300)**3
end where
where(Teff>30000)
  kreac=-1.74046e-11_wp &
    +3.02328e-12_wp*100 &
    -2.39214e-15_wp*100**2 &
    -4.02394e-17_wp*100**3
end where

betanow=kreac*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + NO --> NO+ + O
Teff=30/(16+30._wp)*(16*amu/3/kB*(dv2) &
       + Ts(1:lx1,1:lx2,1:lx3,1) -Tn) + Tn

where (Teff<=3800)
  kreac=6.40408e-13_wp &
    -1.33888e-13_wp*(Teff/300) &
    +7.65103e-14_wp*(Teff/300)**2 &
    -3.11509e-15_wp*(Teff/300)**3 &
    +6.62374e-17_wp*(Teff/300)**4
end where
where (Teff>3800 .and. Teff<=30000)
  kreac=-7.48312e-13_wp &
    +2.31502e-13_wp*(Teff/300) &
    +3.07160e-14_wp*(Teff/300)**2 &
    -2.65436e-16_wp*(Teff/300)**3 &
    +7.76665e-19_wp*(Teff/300)**4
end where
where (Teff>30000)
  kreac=-7.48312e-13_wp &
    +2.31502e-13_wp*(100) &
    +3.07160e-14_wp*(100)**2 &
    -2.65436e-16_wp*(100)**3 &
    +7.76665e-19_wp*(100)**4
end where

betanow=kreac*nn(:,:,:,6)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!O+ + e --> O + hv
betanow=3.7e-12_wp*(250/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.7*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!N2+ + O --> O+ + N2
Teff=16/(28+16._wp)*(28*amu/3/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) -Tn) + Tn

where (Teff <= 1500)
  kreac=1e-11_wp*(300/Teff)**0.23
elsewhere
  kreac=3.6e-12_wp*(300/Teff)**(-0.41)
end where

betanow=kreac*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N+ + O --> O+ + N
betanow=5e-13_wp*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!H+ + O --> O+ + H
Teff=Ts(1:lx1,1:lx2,1:lx3,6)
betanow = (6.e-10_wp)*(8/9._wp)*(((Teff+Tn/4)/(Tn+Teff/16))**0.5)*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,1)=Pr(:,:,:,1)+betanow*ns(1:lx1,1:lx2,1:lx3,6)
Lo(:,:,:,6)=Lo(:,:,:,6)+betanow


!O+ + H --> H+ + O
betanow = 6.0e-10_wp*nn(:,:,:,4)*1e-6_wp
Pr(:,:,:,6)=Pr(:,:,:,6)+betanow*ns(1:lx1,1:lx2,1:lx3,1)
Lo(:,:,:,1)=Lo(:,:,:,1)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! NO+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O+ + NO --> NO+ + O Above


!O2+ + N2 --> NO+ + NO
betanow=5e-16_wp*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!O2+ + N --> NO+ + O
betanow=1.2e-10_wp*nn(:,:,:,5)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!O2+ + NO --> NO+ + O2
betanow=4.6e-10_wp*nn(:,:,:,6)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,4)
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!N2+ + O --> NO+ + N
Teff=16/(28+16._wp)*(28*amu/3/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) - Tn) + Tn

where (Teff <= 1500)
  kreac=1.4e-10_wp*(300/Teff)**0.44
elsewhere
  kreac=5.2e-11_wp*(300/Teff)**(-0.2)
end where

betanow=kreac*nn(:,:,:,1)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N2+ + NO --> NO+ + N2
betanow=4.1e-10_wp*nn(:,:,:,6)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N+ + O2 --> NO+ + O
betanow=2.6e-10_wp*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,2)=Pr(:,:,:,2)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!NO+ + e --> N + O
betanow=4.2e-7_wp*(300/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.85*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,2)=Lo(:,:,:,2)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! N2+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!N2+ + O2 --> O2+ + N2
Teff=32/(28+32._wp)*(28*amu/3/kB*(dv2) &
       +Ts(1:lx1,1:lx2,1:lx3,3) - Tn) + Tn

betanow=5e-11_wp*(300/Teff)*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,3)
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!N2+ + O --> NO+ + N   Above


!N2+ + O --> O+ + N2    Above


!N2+ + O --> NO+ + N    Above


!N2+ + e --> N + N
betanow=1.8e-7_wp*(300/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.39*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,3)=Lo(:,:,:,3)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! O2+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!O2+ + NO --> NO+ + O2  Above


!O+ + O2 --> O2+ + O    Above


!N2+ + O2 --> O2+ + N2  Above


!N+ + O2 --> O2+ + N
betanow=3.1e-10_wp*nn(:,:,:,3)*1e-6_wp
Pr(:,:,:,4)=Pr(:,:,:,4)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!O2+ + e- --> O + O
where (Te <= 1200)
  kreac=1.95e-7_wp* (300/Te)**0.70! See idl code. this may need another te term
elsewhere
  kreac=7.38e-8_wp*(1200/Te)**0.56! See idl code. this may need another te term
end where

betanow=kreac*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,4)=Lo(:,:,:,4)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! N+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!N+ + O --> O+ + N  Above


!N+ + O2 --> NO+ + O Above


!N+ + O2 --> O2+ + N    Above


!N+ + H --> H+ + N
betanow = 3.6e-12_wp*nn(:,:,:,4)*1e-6_wp
Pr(:,:,:,6)=Pr(:,:,:,6)+betanow*ns(1:lx1,1:lx2,1:lx3,5)
Lo(:,:,:,5)=Lo(:,:,:,5)+betanow


!!!!!!!!!!!!!!!!!!!!!!!!!!! H+ REACTIONS !!!!!!!!!!!!!!!!!!!!!!
!H+ + O --> O+ + H above


!O+ + H --> H+ + O above


!N+ + H --> H+ + N above


!H+ + e --> H + hv
betanow=3.7e-12_wp*(250/Ts(1:lx1,1:lx2,1:lx3,lsp))**0.7*ns(1:lx1,1:lx2,1:lx3,lsp)*1e-6_wp
Lo(:,:,:,6)=Lo(:,:,:,6)+betanow

end subroutine srcsContinuity


subroutine srcsMomentum_curv(nn,vn1,Tn,ns,vs1,vs2,vs3,Ts,E1,Q,x,Pr,Lo)

!------------------------------------------------------------
!-------POPULATE SOURCE/LOSS ARRAYS FOR MOMENTUM EQUATION.  ION
!-------PARAMETER ARGUMENTS (AND GRID STUFF) SHOULD INCLUDE GHOST CELLS
!-------NOTE THAT THIS IS THE ONLY SOURCE SUBPROGRAM WHOSE CODE
!-------DIFFERS FROM CARTESIAN TO CURVILINEAR DUE TO PRESSURE
!-------GRADIENT.
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: vn1,Tn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts
real(wp), dimension(-1:,-1:,-1:), intent(in) :: E1
real(wp), dimension(:,:,:,:), intent(in) :: Q
class(curvmesh), intent(in) :: x

real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(inout) :: Pr,Lo
!! intent(out)

integer :: lx1,lx2,lx3,isp,isp2
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: nu,Phisj,Psisj
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: pressure,gradlp1,Epol1,gradQ
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: h1h2h3
real(wp), dimension(0:size(Ts,1)-3,size(Ts,2)-4,size(Ts,3)-4) :: tmpderiv
real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: dh2dx1,dh3dx1,geom
real(wp), dimension(size(E1,1)-4,size(E1,2)-4,size(E1,3)-4) :: E1filt

integer :: ix1,ix2,ix3


lx1=size(Ts,1)-4
lx2=size(Ts,2)-4
lx3=size(Ts,3)-4

Pr=0
Lo=0


!CALCULATE COMMON GEOMETRIC FACTORS USED IN EACH OF THE SPECIES CALCULATIONS
h1h2h3=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
tmpderiv=grad3D1(x%h2(0:lx1+1,1:lx2,1:lx3),x,0,lx1+1,1,lx2,1,lx3)
dh3dx1=tmpderiv(1:lx1,1:lx2,1:lx3)
tmpderiv=grad3D1(x%h3(0:lx1+1,1:lx2,1:lx3),x,0,lx1+1,1,lx2,1,lx3)
dh2dx1=tmpderiv(1:lx1,1:lx2,1:lx3)


!AMBIPOLAR ELECTRIC FIELD
pressure=ns(1:lx1,1:lx2,1:lx3,lsp)*kB*Ts(1:lx1,1:lx2,1:lx3,lsp)
gradlp1=grad3D1(log(pressure),x,1,lx1,1,lx2,1,lx3)
Epol1=kB*Ts(1:lx1,1:lx2,1:lx3,lsp)/qs(lsp)*gradlp1


!THE FIELD INTEGRATED SOLVE ELECTRIC FIELDS ARE NOT RELIABLE BELOW 100KM - AT LEAST NOT ENOUGH TO USE IN THIS CALCULATION
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if (x%alt(ix1,ix2,ix3)<100e3_wp) then
        E1filt(ix1,ix2,ix3)=0
      else
        E1filt(ix1,ix2,ix3)=E1(ix1,ix2,ix3)
      end if
    end do
  end do
end do


do isp=1,lsp
  !ION-NEUTRAL COLLISIONS
  do isp2=1,ln
    call maxwell_colln(isp,isp2,nn,Tn,Ts,nu)

    Lo(:,:,:,isp)=Lo(:,:,:,isp)+nu
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*nu*vn1
  end do


  !ION-ION
  do isp2=1,lsp
    call coulomb_colln(isp,isp2,ns,Ts,vs1,nu,Phisj,Psisj)

    Lo(:,:,:,isp)=Lo(:,:,:,isp)+nu*Phisj
    Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp) &
                  *nu*Phisj*vs1(1:lx1,1:lx2,1:lx3,isp2)
  end do


  !ION PRESSURE
  pressure=ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)
  gradlp1=grad3D1(log(pressure),x,1,lx1,1,lx2,1,lx3)                         !derivative should be from 1:lx1
  !might need to limit the gradient to non-null points like 2D MATLAB code


  !ARTIFICIAL VISCOSITY
  gradQ=grad3D1(Q(:,:,:,isp),x,1,lx1,1,lx2,1,lx3)                         !derivative should be from 1:lx1

  !GEOMETRIC FACTORS ARISING FROM ADVECTINO OF 1-COMPONENT OF MOMENTUM DENSITY
  geom=(vs2(1:lx1,1:lx2,1:lx3,isp)**2*x%h3(1:lx1,1:lx2,1:lx3)*dh2dx1+ &
        vs3(1:lx1,1:lx2,1:lx3,isp)**2*x%h2(1:lx1,1:lx2,1:lx3)*dh3dx1)*ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/h1h2h3


  !ACCUMULATED ALL FORCES
!      Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)*(E1+Epol1) &
  Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)*(E1filt+Epol1) &
                -pressure*gradlp1 &
                -gradQ &
                +geom &
                +ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*x%g1

  ! print*, minval(E1filt),maxval(E1filt),minval(E1),maxval(E1)
!  Pr(:,:,:,isp)=ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)*(E1filt+Epol1) &
!                -pressure*gradlp1 &
!                -gradQ &
!                +geom &
!                +ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*x%g1
end do

end subroutine srcsMomentum_curv


subroutine srcsEnergy(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo,E2,E3,x,cfg)
  !------------------------------------------------------------
  !-------POPULATE SOURCE/LOSS ARRAYS FOR ENERGY EQUATION.  ION
  !-------PARAMETER ARGUMENTS SHOULD INCLUDE GHOST CELLS
  !------------------------------------------------------------
  
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,vs1,vs2,vs3,Ts
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: E2,E3
  class(curvmesh), intent(in) :: x !Added for FBI, need BMAG
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp), intent(inout) :: Pr,Lo
  !! intent(out)
  type(gemini_cfg), intent(in) :: cfg
  integer :: lx1,lx2,lx3,isp,isp2
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: nu,Phisj,Psisj
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: fact,iePT,ieLT,f,g    !work arrays
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: N2vibrationalLoss, O2vibrationalLoss
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: FBIproduction,FBIlossfactor    !FBI array
  real(wp) :: sfact
 

  lx1=size(Ts,1)-4
  lx2=size(Ts,2)-4
  lx3=size(Ts,3)-4
  
  Pr=0
  Lo=0
  iePT=0
  ieLT=0
  
  
  !ELASTIC COLLISIONS
  do isp=1,lsp
    !ION-NEUTRAL
    do isp2=1,ln
      call maxwell_colln(isp,isp2,nn,Tn,Ts,nu)
  
      !HEAT TRANSFER
      fact=2*nu/(ms(isp)+mn(isp2))
      Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*kB/(gammas(isp)-1)*fact*Tn
      Lo(:,:,:,isp)=Lo(:,:,:,isp)+ms(isp)*fact
  
  
      !FRICTION
      fact=fact*mn(isp2)/3
      Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/(gammas(isp)-1) &
                    *((vs1(1:lx1,1:lx2,1:lx3,isp)-vn1)**2+(vs2(1:lx1,1:lx2,1:lx3,isp)-vn2)**2 &
                    +(vs3(1:lx1,1:lx2,1:lx3,isp)-vn3)**2)*fact     !vn's should be correct shape for this...
    end do
  
    !ION-ION
    do isp2=1,lsp
      call coulomb_colln(isp,isp2,ns,Ts,vs1,nu,Phisj,Psisj)
  
      !HEAT TRANSFER
      fact=2*nu*Psisj/(ms(isp)+ms(isp2))
      Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*kB/(gammas(isp)-1) &
                    *fact*Ts(1:lx1,1:lx2,1:lx3,isp2)
      Lo(:,:,:,isp)=Lo(:,:,:,isp)+ms(isp)*fact
  
      !FRICTION
  !        fact=2*nu*Phisj/(ms(isp)+ms(isp2))*mn(isp2)/3     !this is the error that was causing the runtime problem with -O3 on phys_consts.f90.  Much thanks to Guy Grubbs for finding this longstanding error.
      fact=2*nu*Phisj/(ms(isp)+ms(isp2))*ms(isp2)/3
      Pr(:,:,:,isp)=Pr(:,:,:,isp)+ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)/(gammas(isp)-1) &
                    *((vs1(1:lx1,1:lx2,1:lx3,isp)-vs1(1:lx1,1:lx2,1:lx3,isp2))**2 &
                     +(vs2(1:lx1,1:lx2,1:lx3,isp)-vs2(1:lx1,1:lx2,1:lx3,isp2))**2 &
                     +(vs3(1:lx1,1:lx2,1:lx3,isp)-vs3(1:lx1,1:lx2,1:lx3,isp2))**2)*fact
    end do
  end do
  
  
  !INELASTIC COLLISIONS FOR ELECTRONS, ROTATIONAL
  if (cfg%flagevibcool>0) then
    sfact=elchrg/kB*(gammas(lsp)-1);   !cf. S&N 2010, electron energy equatoin section fixed by JMDP
    nu=sfact*5.2e-15_wp*nn(:,:,:,3)*1e-6_wp/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !O2 rotational excitation 5.2e-15
    iePT=nu*Tn
    ieLT=nu
    nu=sfact*3.5e-14_wp*nn(:,:,:,2)*1e-6_wp/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !N2 rot. exc. 3.5e-14
    iePT=iePT+nu*Tn
    ieLT=ieLT+nu
    
    call N2vib(nn,Tn,Ts,N2vibrationalLoss)
    call O2vib(nn, Tn, Ts, O2vibrationalLoss)
    iePT=iePT-max(O2vibrationalLoss,0._wp)-max(N2vibrationalLoss,0._wp)
  else
    !INELASTIC COLLISIONS FOR ELECTRONS, ROTATIONAL
    sfact=elchrg/kB*(gammas(lsp)-1);   !cf. S&N 2010, electron energy equatoin section
    nu=sfact*6.9e-14_wp*nn(:,:,:,3)*1e-6_wp/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !O2 rotational excitation
    iePT=nu*Tn
    ieLT=nu
    nu=sfact*2.9e-14_wp*nn(:,:,:,2)*1e-6_wp/sqrt(Ts(1:lx1,1:lx2,1:lx3,lsp))    !N2 rot. exc.
    iePT=iePT+nu*Tn;
    ieLT=ieLT+nu;
     
    f=1.06e4_wp+7.51e3_wp*tanh(1.10e-3_wp*(Ts(1:lx1,1:lx2,1:lx3,lsp)-1800))
    g=3300+1.233_wp*(Ts(1:lx1,1:lx2,1:lx3,lsp)-1000)-2.056e-4_wp &
      *(Ts(1:lx1,1:lx2,1:lx3,lsp)-1000)*(Ts(1:lx1,1:lx2,1:lx3,lsp)-4000)
    fact=sfact*2.99e-12_wp*nn(:,:,:,2)*1e-6_wp*exp(f*(Ts(1:lx1,1:lx2,1:lx3,lsp)-2000) &
            /Ts(1:lx1,1:lx2,1:lx3,lsp)/2000)*(exp(-g*(Ts(1:lx1,1:lx2,1:lx3,lsp)-Tn) &
            /Ts(1:lx1,1:lx2,1:lx3,lsp)/Tn)-1)    !N2 vibrational excitation
    iePT=iePT-max(fact,0._wp);
     f=3300-839*sin(1.91e-4_wp*(Ts(1:lx1,1:lx2,1:lx3,lsp)-2700))
     fact=sfact*5.196e-13_wp*nn(:,:,:,3)*1e-6_wp*exp(f*(Ts(1:lx1,1:lx2,1:lx3,lsp)-700) &
      /Ts(1:lx1,1:lx2,1:lx3,lsp)/700)*(exp(-2770*(Ts(1:lx1,1:lx2,1:lx3,lsp)-Tn) &
      /Ts(1:lx1,1:lx2,1:lx3,lsp)/Tn)-1)    !O2 vibrational excitation
    iePT=iePT-max(fact,0._wp);    
  end if
 
  
  !This includes losses of the FBI part
  !CORRECT TEMP EXPRESSIONS TO CORRESPOND TO INTERNAL ENERGY SOURCES
  !This would be the place to include FBI heating probably just add to iePT
  if (cfg%flagFBI>0) then 
    call FBIheating(nn,Tn,ns,Ts,E2,E3,x,FBIproduction,FBIlossfactor)
    Pr(:,:,:,lsp)=Pr(:,:,:,lsp)+FBIproduction+(iePT*FBIlossfactor)*ns(1:lx1,1:lx2,1:lx3,lsp)*kB/(gammas(lsp)-1)   !Arg, forgot about the damn ghost cells in original code...
    Lo(:,:,:,lsp)=Lo(:,:,:,lsp)+(ieLT*FBIlossfactor)
  else  
  ! !This is for no FBI simulation
    Pr(:,:,:,lsp)=Pr(:,:,:,lsp)+iePT*ns(1:lx1,1:lx2,1:lx3,lsp)*kB/(gammas(lsp)-1)   !Arg, forgot about the damn ghost cells in original code...
    Lo(:,:,:,lsp)=Lo(:,:,:,lsp)+ieLT
  end if
end subroutine srcsEnergy



subroutine O2vib(nn,Tn,Ts,O2VibrationalLoss)
  real(wp), dimension(:,:,:,:), intent(in) :: nn !Neutral density
  real(wp), dimension(:,:,:), intent(in) :: Tn !neutral temperature
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: Ts !Plasma density and temperature
  real(wp), dimension(:,:,:), intent(out) :: O2VibrationalLoss
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) ::LogQTe, QTe, Te, Teaux
  integer :: lx1,lx2,lx3
  
  lx1=size(Ts,1)-4
  lx2=size(Ts,2)-4
  lx3=size(Ts,3)-4
  
  !Define Te, no ghost cells
  Te=Ts(1:lx1,1:lx2,1:lx3,lsp)
  Teaux=3800.0_wp
  
  !! Calculate Log10(Q(Te))
  where (Te<=Teaux)
  LogQTe = -19.9171_wp &
           +0.0267_wp*Te &
           -3.9960e-5_wp*Te**2 &
           +3.5187e-8_wp*Te**3 &
           -1.9228e-11_wp*Te**4 &
           +6.6865e-15_wp*Te**5 &
           -1.4791e-18_wp*Te**6 &
           +2.0127e-22_wp*Te**7 &
           -1.5346e-26_wp*Te**8 &
           +5.0148e-31_wp*Te**9 
  elsewhere
  LogQTe = -19.9171_wp &
           +0.0267_wp*Teaux &
           -3.9960e-5_wp*Teaux**2 &
           +3.5187e-8_wp*Teaux**3 &
           -1.9228e-11_wp*Teaux**4 &
           +6.6865e-15_wp*Teaux**5 &
           -1.4791e-18_wp*Teaux**6 &
           +2.0127e-22_wp*Teaux**7 &
           -1.5346e-26_wp*Teaux**8 &
           +5.0148e-31_wp*Teaux**9 
  end where
  
  !Because the loss factor is a fitting of the logarithmic base 10 value of it. Multiply by LOG10 to change to natural log
  QTe = EXP(LogQTe*LOG(10.0_wp)) !Make it linear 
  O2VibrationalLoss=nn(:,:,:,3)*1.0e-6_wp*QTe*(1-EXP(2239.0_wp*((Tn-Te)/(Te*Tn))))
  !print*,nn(:,:,:,3)
end subroutine O2vib


subroutine N2vib(nn,Tn,Ts,N2VibrationalLoss)
  real(wp), dimension(:,:,:,:), intent(in) :: nn !Neutral density
  real(wp), dimension(:,:,:), intent(in) :: Tn !neutral temperature
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: Ts !Plasma density and temperature
  
  real(wp), dimension(:,:,:), intent(out) :: N2VibrationalLoss
  
  integer :: ix1,ix2,ix3,lx1,lx2,lx3, isp
  
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: vibloss, Te
  
  real(wp), parameter :: E1 = 3353.0_wp
  real(wp) :: LQ0Te=0.0_wp, Q0Te=0.0_wp, LQ1Te=0.0_wp, Q1Te=0.0_wp, &
              Tei=0.0_wp, Tni=0.0_wp, nni=0.0_wp, &
              STerm0, STerm1
  
  real(wp), parameter :: A0(10)=[real(wp) :: -2.025_wp, &
                                              7.066_wp, & 
                                              8.211_wp, &
                                              9.713_wp, &
                                              10.353_wp, &
                                              10.819_wp, &
                                              10.183_wp, &
                                              12.698_wp, &
                                              14.710_wp, &
                                              17.538_wp]*(-1.0_wp)
  
  real(wp), parameter :: B0(10)=[real(wp) :: 8.782e-2_wp, &
                                             1.001_wp, &
                                             1.092_wp, &
                                             1.204_wp, &
                                             1.243_wp, &
                                             1.244_wp, &
                                             1.185_wp, &
                                             1.309_wp, &
                                             1.409_wp, &
                                             1.600_wp]*(1.0e-2_wp)
  
  real(wp), parameter :: C0(10)=[real(wp) :: -2.954e-1_wp, &
                                              3.066_wp, &
                                              3.369_wp, &
                                              3.732_wp, &
                                              3.850_wp, &
                                              3.771_wp, &
                                              3.570_wp, &
                                              3.952_wp, &
                                              4.249_wp, &
                                              4.916_wp]*(-1.0e-6_wp)
  
  real(wp), parameter :: D0(10)=[real(wp) :: -9.562e-1_wp, &
                                              4.436_wp, &
                                              4.891_wp, &
                                              5.431_wp, &
                                              5.600_wp, &
                                              5.385_wp, &
                                              5.086_wp, &
                                              5.636_wp, &
                                              6.058_wp, &
                                              7.128_wp]*(1.0e-10_wp)
  
  real(wp), parameter :: F0(10)=[real(wp) :: -7.252e-1_wp, &
                                              2.449_wp, &
                                              2.706_wp, &
                                              3.008_wp, &
                                              3.100_wp, &
                                              2.936_wp, &
                                              2.769_wp, &
                                              3.071_wp, &
                                              3.300_wp, &
                                              3.941_wp]*(-1.0e-14_wp)
  
  real(wp), parameter :: A1(8)=[real(wp) :: 3.413_wp, &
                                             4.160_wp, &
                                             5.193_wp, &
                                             5.939_wp, &
                                             8.261_wp, &
                                             8.185_wp, &
                                             10.823_wp, &
                                             11.273_wp]*(-1.0_wp)
  
  real(wp), parameter :: B1(8)=[real(wp) :: 7.326e-1_wp, &
                                            7.803e-1_wp, &
                                            8.360e-1_wp, &
                                            8.807e-1_wp, &
                                            1.010_wp, &
                                            1.010_wp, &
                                            1.199_wp, &
                                            1.283_wp]*(1.0e-2_wp)
  
  real(wp), parameter :: C1(8)=[real(wp) :: 2.200_wp, &
                                            2.352_wp, &
                                            2.526_wp, &
                                            2.669_wp, &
                                            3.039_wp, &
                                            3.039_wp, &
                                            3.620_wp, &
                                            3.879_wp]*(-1.0e-6_wp)
  
  real(wp), parameter :: D1(8)=[real(wp) :: 3.128_wp, &
                                            3.352_wp, &
                                            3.606_wp, &
                                            3.806_wp, &
                                            4.318_wp, &
                                            4.318_wp, &
                                            5.159_wp, &
                                            5.534_wp]*(1.0e-10_wp)
  
  real(wp), parameter :: F1(8)=[real(wp) :: 1.702_wp, &
                                            1.828_wp, &
                                            1.968_wp, &
                                            2.073_wp, &
                                            2.347_wp, &
                                            2.347_wp, &
                                            2.810_wp, &
                                            3.016_wp]*(-1.0e-14_wp)
  
                                            
  
  lx1=size(Ts,1)-4
  lx2=size(Ts,2)-4
  lx3=size(Ts,3)-4
  
  !Define Te, no ghost cells
  Te=Ts(1:lx1,1:lx2,1:lx3,lsp)
  
  !loop over all indexes in the grid
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
  
        STerm0=0.0_wp !make summation terms 0 again, previous cell will have numbers here
        STerm1=0.0_wp
        Tei=Te(ix1,ix2,ix3) !Define it for selection
        Tni=Tn(ix1,ix2,ix3)
        nni=nn(ix1,ix2,ix3,2) !just N2
  
        if (Tei<=1500) then       !Case 1, lower temperature
          LQ0Te=-6.462_wp+3.151e-2_wp*Tei-4.075e-5_wp*Tei**2+2.439e-8_wp*Tei**3-5.479e-12_wp*Tei**4-16.0_wp;
          Q0Te=EXP(LQ0Te*LOG(10.0_wp))
          vibloss(ix1,ix2,ix3)=nni*1e-6_wp*(1-EXP(-E1/Tni))*Q0Te*(1-exp(E1*((Tni-Tei)/(Tei*Tni))))
        else if (Tei<=6000) then !Case 2, have to calculate summations
          do isp=1,10
            LQ0Te=A0(isp)+Tei*B0(isp)+Tei**2*C0(isp)+Tei**3*D0(isp)+Tei**4*F0(isp)-16.0_wp;
            Q0Te=EXP(LQ0Te*LOG(10.0_wp))
            STerm0=STerm0+Q0Te*(1-exp(isp*E1*((Tni-Tei)/(Tei*Tni))))
          end do
          do isp=1,8
            LQ1Te=A1(isp)+Tei*B1(isp)+Tei**2*C1(isp)+Tei**3*D1(isp)+Tei**4*F1(isp)-16.0_wp;
            Q1Te=EXP(LQ1Te*LOG(10.0_wp))
            STerm1=STerm1+Q1Te*(1-exp(isp*E1*((Tni-Tei)/(Tei*Tni))))
          end do
          vibloss(ix1,ix2,ix3)=nni*1e-6_wp*(1-EXP(-E1/Tni))*Sterm0+ &
                               nni*1e-6_wp*(1-EXP(-E1/Tni))*EXP(-E1/Tni)*Sterm0
        else 
          Tei=6000 !the to max temperature value and do the same as above
          do isp=1,10
            LQ0Te=A0(isp)+Tei*B0(isp)+Tei**2*C0(isp)+Tei**3*D0(isp)+Tei**4*F0(isp)-16.0_wp;
            Q0Te=EXP(LQ0Te*LOG(10.0_wp))
            STerm0=STerm0+Q0Te*(1-exp(isp*E1*((Tni-Tei)/(Tei*Tni))))
          end do
          do isp=1,8
            LQ1Te=A1(isp)+Tei*B1(isp)+Tei**2*C1(isp)+Tei**3*D1(isp)+Tei**4*F1(isp)-16.0_wp;
            Q1Te=EXP(LQ1Te*LOG(10.0_wp))
            STerm1=STerm1+Q1Te*(1-exp(isp*E1*((Tni-Tei)/(Tei*Tni))))
          end do
          vibloss(ix1,ix2,ix3)=nni*1e-6_wp*(1-EXP(-E1/Tni))*Sterm0+ &
                               nni*1e-6_wp*(1-EXP(-E1/Tni))*EXP(-E1/Tni)*Sterm0
        end if   
      end do
    end do
  end do
  
  N2VibrationalLoss=vibloss
end subroutine N2vib



subroutine FBIheating(nn,Tn,ns,Ts,E2,E3,x,FBIproduction,FBIlossfactor)
  !! Inputs Needed
  real(wp), dimension(:,:,:,:), intent(in) :: nn !Neutral density
  real(wp), dimension(:,:,:), intent(in) :: Tn !neutral temperature
  real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts !Plasma density and temperature
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: E2,E3 !Electric Field
  class(curvmesh), intent(in) :: x !Grid, doing this because BMAG is stored here
  
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: FBIproduction,FBIlossfactor ! Two terms, one is heating and the other one is a factor for cooling. 
  
  !!Internal Arrays
  integer :: isp,isp2,lx1,lx2,lx3
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp) :: nsuAvg 
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,lsp-1) :: niW 
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,ln) :: nuW 
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4,2) :: nuAvg, msAvg, TsAvg 
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: Bmagnitude, nu, nsAvg, omegae, omegai, ki, ke, phi
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: Eth0, Ethreshold, Emagnitude
  real(wp), dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: heatingfirst, heatingsecond, heatingtotal, lossfactor
  integer, dimension(size(Ts,1)-4,size(Ts,2)-4,size(Ts,3)-4) :: FBIbinary
  
  !real(wp), dimension(lbound(Ts,1):ubound(Ts,1),lbound(Ts,2):ubound(Ts,2),lbound(Ts,3):ubound(Ts,3),lsp) :: Tsfix
  !real(wp) :: TMAX
  
  lx1=x%lx1
  lx2=x%lx2
  lx3=x%lx3
  
  !Bmagnitude=x%Bmag(1:lx1,1:lx2,1:lx3)
  Bmagnitude=x%Bmag(1:lx1,1:lx2,1:lx3)
  Emagnitude=sqrt(E2(1:lx1,1:lx2,1:lx3)**2+E3(1:lx1,1:lx2,1:lx3)**2) !!Already evaluated with no ghost cells
  
  !!Initialize arrays a 0s and loss as 1s
  nuAvg=0.0_wp
  nsuAvg=0.0_wp
  msAvg=0.0_wp
  nsAvg=0.0_wp
  TsAvg=0.0_wp
  FBIproduction=0.0_wp
  FBIlossfactor=0.0_wp
  lossfactor=1.0_wp
  heatingfirst=0.0_wp
  heatingsecond=0.0_wp
  heatingtotal=0.0_wp
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
  !nsAvg=ns(1:lx1,1:lx2,1:lx3,lsp) !! Assume qneutrality, could be wrong
  !! Average just O2+ and NO+
  nsAvg=(ns(1:lx1,1:lx2,1:lx3,2)*niW(:,:,:,2)+ns(1:lx1,1:lx2,1:lx3,4)*niW(:,:,:,4))/(niW(:,:,:,2)+niW(:,:,:,4))
  
  !! ki value
  omegai=elchrg*Bmagnitude/msAvg(:,:,:,1) !! Would this work?, it will, I defined Bmagnitude above
  ki=abs(omegai/nuAvg(:,:,:,1)) !! Could do ABS to be sure of the sign
  !! ke value
  omegae=elchrg*Bmagnitude/msAvg(:,:,:,2) 
  ke=abs(omegae/nuAvg(:,:,:,2))
  
  !!Phi value 1/(ki*ke)
  phi=1.0_wp/(ke*ki)
  
  !!Average ion temperature
  TsAvg(:,:,:,1)=(Ts(1:lx1,1:lx2,1:lx3,2)*niW(:,:,:,2)+Ts(1:lx1,1:lx2,1:lx3,4)*niW(:,:,:,4))/(niW(:,:,:,2)+niW(:,:,:,4))
  TsAvg(:,:,:,2)=Ts(1:lx1,1:lx2,1:lx3,lsp)
  
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
  
  !Calculate heating term only where FBI is possible
  where (FBIbinary==1)
    heatingfirst=(msAvg(:,:,:,1)*nuAvg(:,:,:,1)*nsAvg*(ki**2)*(Emagnitude-Ethreshold)**2)/((1.0_wp+ki**2)*Bmagnitude**2)
    heatingsecond=((Emagnitude/Ethreshold)*(1.0_wp+phi)-1.0_wp)
    heatingtotal=heatingfirst*heatingsecond
  end where
  
  call LossFactorCalc(TsAvg(:,:,:,2),FBIbinary,lossfactor)
  
  FBIproduction=heatingtotal
  FBIlossfactor=lossfactor
end subroutine FBIheating


subroutine LossFactorCalc(Te,FBIBinary,FourierLossFactor)
  !! Inputs Needed
  real(wp), dimension(:,:,:), intent(in) :: Te!Temperature, without ghost cells
  integer, dimension(:,:,:), intent(in) :: FBIbinary
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: FourierLossFactor !  
  !!Internal Arrays
  real(wp), dimension(size(Te,1),size(Te,2),size(Te,3)) :: costerms, sinterms, a0
  
  real(wp), parameter ::  a1 = -159.3814_wp, &      !(-177.7, -141.1)
                          a2 = 1.7150e+03_wp, &     !(1097, 2333)
                          a3 = 244.4458_wp, &       !(216.4, 272.5)
                          a4 = -579.3561_wp, &       !(-796.1, -362.6)
                          a5 = -94.0672_wp, &       !(-104.9, -83.24)
                          a6 = 75.7883_wp, &        !(45.45, 106.1)
                          a7 = 9.7414_wp, &         !(8.655, 10.83)
                          a8 = -1.7097_wp        !(-2.494, -0.9256)
  
  real(wp), parameter ::  b1 = 2.2211e+03_wp, &     !%(1427, 3015)
                          b2 = 249.2960_wp, &       !%(220.7, 277.9)
                          b3 = -1.1022e+03_wp, &    !%(-1506, -698.8)
                          b4 = -175.0776_wp, &      !%(-195.2, -155)
                          b5 = 241.7140_wp, &       !%(148.6, 334.8)
                          b6 = 36.9935_wp, &        !%(32.77, 41.22)
                          b7 = -15.9934_wp, &       !%(-22.75, -9.238)
                          b8 = -1.3170_wp        !%(-1.457, -1.177)
  
  real(wp), parameter ::  w = 1.0402e-04_wp         !%(0.0001018, 0.0001062)
  
  a0 = -1.2115e+03_wp
  costerms=0.0_wp
  sinterms=0.0_wp
  
  where (FBIBinary==1)
    costerms = a1*cos(Te*w) + a2*cos(2.0_wp*Te*w) + a3*cos(3.0_wp*Te*w) + a4*cos(4.0_wp*Te*w) &
              + a5*cos(5.0_wp*Te*w) + a6*cos(6*Te*w) + a7*cos(7.0_wp*Te*w) + a8*cos(8.0_wp*Te*w)
    sinterms = b1*sin(Te*w) + b2*sin(2.0_wp*Te*w) + b3*sin(3.0_wp*Te*w) + b4*sin(4.0_wp*Te*w) &
              + b5*sin(5.0_wp*Te*w) + b6*sin(6*Te*w) + b7*sin(7.0_wp*Te*w) + b8*sin(8.0_wp*Te*w)
    !Because the loss factor is a fitting of the logarithmic base 10 value of it. Multiply by LOG10 to change to natural log
    FourierLossFactor = EXP((a0 + costerms + sinterms)*LOG(10.0_wp)) !Make it linear 
  elsewhere 
    FourierLossFactor = 1.0_wp
  end where
  
  where (Te>=30000.0_wp)
    FourierLossFactor=0.003811931223844_wp
  end where
  
  where (Te<=600.0_wp)
    FourierLossFactor=1.0_wp
  end where
end subroutine LossFactorCalc



end module sources
