module diffusion

!! This module sets up the ionospheric diffusion problem and then passes it off to the parabolic solvers.

use phys_consts, only: kB,lsp,gammas,mindensdiv, wp
!! sizes are not imported in case we want to do subgrid diffusion
use meshobj, only : curvmesh
use grid, only : gridflag
use PDEparabolic, only : backEuler1D,TRBDF21D

implicit none (type, external)
private
public :: diffusion_prep, TRBDF23D, backEuler3D

interface TRBDF23D
  module procedure TRBDF23D_curv
end interface TRBDF23D

interface backEuler3D
  module procedure backEuler3D_curv
end interface backEuler3D

integer, dimension(2), protected :: BCtype=[0,0]

contains
  impure subroutine diffusion_prep(isp,x,lambda,betacoeff,ns,T,A,B,C,D,E,Tn,Teinf)
    !! COMPUTE COEFFICIENTS IN DIFFUSION EQUATIONS AND LOAD UP GHOST CELLS
    !!
    !! Note: done on a per species basis. This is never called over the full grid
    integer, intent(in) :: isp
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:), intent(in) :: lambda,betacoeff
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: ns
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: T
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), &
             intent(out) :: A,B,C,D,E
    real(wp), dimension(:,:,:), intent(in) :: Tn
    real(wp), intent(in) :: Teinf
    real(wp) :: Tn0
    integer :: lx1,lx2,lx3,ix2,ix3
    
    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4
    
    !COEFFICIENTS OF PARABOLIC EQUAITON
    A(:,:,:)=0._wp
     C(:,:,:)=(gammas(isp)-1._wp)/kB/max(ns(1:lx1,1:lx2,1:lx3),mindensdiv)/ &
               (x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3))
    B(:,:,:)=C(:,:,:)*betacoeff/x%h1(1:lx1,1:lx2,1:lx3)    !beta must be set to zero if not electrons!
    D=lambda*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)
    E=0._wp

    ! Determine what type of boundary conditions we are trying to use
    call set_BCtype(Teinf,gridflag)
    
    !SET THE BOUNDARY CONDITIONS BASED ON GRID TYPE
    ! if Neumann need to scale heat flux by thermal conductivity and metric factor...
    if (gridflag==0) then    !closed dipole grid, both ends are thermalized against neutrals
      do ix3=1,lx3
        do ix2=1,lx2
          Tn0=Tn(1,ix2,ix3)
          T(0,ix2,ix3)=Tn0
          Tn0=Tn(lx1,ix2,ix3)
          T(lx1+1,ix2,ix3)=Tn0
        end do
      end do
    else if (gridflag==1) then    !inverted grid, bottom altitude is thermalized to neutrals, top to electrons (possibly heat flow)
      do ix3=1,lx3
        do ix2=1,lx2
          Tn0=Tn(lx1,ix2,ix3)
          T(lx1+1,ix2,ix3)=Tn0      !bottom
          if (BCtype(1)==0) then    ! user wants Dirichlet on 'top'
            T(0,ix2,ix3)=Teinf
          else                      ! Neumann
            if (isp==7) then
              T(0,ix2,ix3)=-Teinf*x%h1(1,ix2,ix3)*x%dx1(2)/lambda(1,ix2,ix3)
            else
              T(0,ix2,ix3)=0._wp         
            end if
          end if
        end do
      end do
    else                          !non-inverted, standard.  Bottom is logical first element of array...
      do ix3=1,lx3
        do ix2=1,lx2
          Tn0=Tn(1,ix2,ix3)
          T(0,ix2,ix3)=Tn0            !bottom
          if (BCtype(2)==0) then      ! Dirichlet on 'top'
            T(lx1+1,ix2,ix3)=Teinf
          else                        ! Neumann
            if (isp==7) then 
              T(lx1+1,ix2,ix3)=-Teinf*x%h1(lx1,ix2,ix3)*x%dx1(lx1)/lambda(lx1,ix2,ix3) !top for neumman
            else
              T(lx1+1,ix2,ix3)=0._wp
            end if
          end if
        end do
      end do
    end if
  end subroutine diffusion_prep


  subroutine set_BCtype(Teinf,gridflag)
    real(wp), intent(in) :: Teinf
    integer, intent(in) :: gridflag

    if (gridflag==1) then
      if (Teinf<1.0) then
        BCtype=[1,0]
      else
        BCtype=[0,0]
      end if
    else if (gridflag==2) then
      if (Teinf<1.0) then           
        BCtype=[0,1]
      else
        BCtype=[0,0]
      end if
    else
      BCtype=[0,0]
    end if
  end subroutine set_BCtype
  
  
  function backEuler3D_curv(f,A,B,C,D,E,dt,x)
    !! SOLVE A 3D SEQUENCE OF 1D DIFFUSION PROBLEMS.
    !! GHOST CELLS ARE ACCOMMODATED AS THEY PROVIDE
    !! A CONVENIENT MEMORY SPACE FOR BOUNDARY CONDITIONS.
    
    real(wp), dimension(:,:,:), intent(in) :: A,B,C,D,E   !trimmed to grid size
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: f    !expected to include ghosts
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    
    integer :: ix2,ix3,lx1,lx2,lx3
    real(wp),dimension(size(f,1)-4) :: fx1slice
    
    real(wp), dimension(-1:size(f,1)-2,-1:size(f,2)-2,-1:size(f,3)-2) :: backEuler3D_curv
    
    lx1=size(f,1)-4
    lx2=size(f,2)-4
    lx3=size(f,3)-4
    
    do ix3=1,lx3
      do ix2=1,lx2
        fx1slice=f(1:lx1,ix2,ix3)
        fx1slice=backEuler1D(fx1slice,A(:,ix2,ix3), &
                      B(:,ix2,ix3),C(:,ix2,ix3),D(:,ix2,ix3),E(:,ix2,ix3), &
                      f(0,ix2,ix3),f(lx1+1,ix2,ix3),dt,BCtype,x%dx1,x%dx1i)
        !! inner ghost cells include boundary conditions
        backEuler3D_curv(1:lx1,ix2,ix3)=fx1slice
      end do
    end do
  end function backEuler3D_curv
  
  
  function TRBDF23D_curv(f,A,B,C,D,E,dt,x)
    !! SOLVE A 3D SEQUENCE OF 1D DIFFUSION PROBLEMS.
    !! GHOST CELLS ARE ACCOMMODATED AS THEY PROVIDE
    !! A CONVENIENT MEMORY SPACE FOR BOUNDARY CONDITIONS.
    !!
    !! Note that this function also plays the role of abstracting
    !! away the grid structure so that the individual 1D lines are
    !! solved using arrays for differences instead of structure members
    
    !> trimmed to grid size
    real(wp), dimension(:,:,:), intent(in) :: A,B,C,D,E
    
    !> expected to include ghosts
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: f
    
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    
    integer :: ix2,ix3,lx1,lx2,lx3
    real(wp),dimension(size(f,1)-4) :: fx1slice
    
    real(wp), dimension(-1:size(f,1)-2,-1:size(f,2)-2,-1:size(f,3)-2) :: TRBDF23D_curv
    
    lx1=size(f,1)-4
    lx2=size(f,2)-4
    lx3=size(f,3)-4
    
    do ix3=1,lx3
      do ix2=1,lx2
        fx1slice=f(1:lx1,ix2,ix3)
        fx1slice=TRBDF21D(fx1slice,A(:,ix2,ix3), &
                      B(:,ix2,ix3),C(:,ix2,ix3),D(:,ix2,ix3),E(:,ix2,ix3), &
                      f(0,ix2,ix3),f(lx1+1,ix2,ix3),dt,BCtype,x%dx1,x%dx1i)
        !! inner ghost cells include boundary conditions
        TRBDF23D_curv(1:lx1,ix2,ix3)=fx1slice
      end do
    end do
  end function TRBDF23D_curv
end module diffusion
