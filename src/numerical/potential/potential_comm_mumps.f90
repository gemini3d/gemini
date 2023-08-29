module potential_comm

!! THIS MODULE IS MEANT TO WORK WITH THE MUMPS 2D
!! INTEGRATED SOLVER IF THE GRID IS 3D, OR A FIELD-RESOLVED
!! SOLVER IF THE GRID IS 2D (MUMPS CAN'T HANDLE 3D VERY WELL).
!!
!! NOTE THAT ONLY THE CURVILINEAR FUNCTION ARE UP TO DATE.

use, intrinsic :: ieee_arithmetic

use phys_consts, only: wp, pi, lsp, debug, ms, qs, kB
use gemini_work_def, only: gemini_work
use grid, only: gridflag, lx1,lx2,lx3,lx2all,lx3all
use meshobj, only: curvmesh
use efielddataobj, only: efielddata
use collisions, only: conductivities, capacitance, NLConductivity
use calculus, only: div3d, integral3d1, grad3d1, grad3d2, grad3d3, integral3d1_curv_alt
use potentialBCs_mumps, only: potentialbcs2D, potentialbcs2D_fileinput, compute_rootBGEfields
use potential_mumps, only: potential3D_fieldresolved_decimate, &
                            potential2D_fieldresolved, &
                            potential2D_polarization, &
                            potential2D_polarization_periodic, &
                            potential3D_fieldresolved, &
                            potential3D_fieldresolved_truncate
use PDEelliptic, only: elliptic_workers
use mpimod, only: mpi_cfg, tag=>gemini_mpi, &
bcast_send, bcast_recv, gather_recv, gather_send, halo, bcast_send3D_ghost, bcast_recv3D_ghost
use gemini3d_config, only: gemini_cfg

use mpi_f08, only: mpi_send, mpi_recv, mpi_integer, mpi_comm_world, mpi_status_ignore

implicit none (type, external)
private
public :: electrodynamics, halo_pot, potential_sourceterms, pot2perpfield, velocities, get_BGEfields, &
            acc_perpconductioncurrents,acc_perpwindcurrents,acc_perpgravcurrents,acc_pressurecurrents, &
            parallel_currents,polarization_currents,BGfields_boundaries_root,BGfields_boundaries_worker, &
            acc_perpBGconductioncurrents

!! overloading to deal with vestigial cartesian->curvilinear code
interface electrodynamics
  module procedure electrodynamics_curv
end interface electrodynamics

interface potential_root_mpi
  module procedure potential_root_mpi_curv
end interface potential_root_mpi

interface ! potential_worker.f90
  module subroutine potential_workers_mpi(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav, &
                                            muP,muH, &
                                            incap,vs2,vs3,vn2,vn3,cfg,B1,ns,Ts,x, &
                                            flagdirich,E02src,E03src,Vminx1slab,Vmaxx1slab, &
                                            E1,E2,E3,J1,J2,J3)
    integer, intent(in) :: it
    real(wp), intent(in) :: t,dt
    real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
    real(wp), dimension(:,:,:), intent(in) ::  incap
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3
    real(wp), dimension(:,:,:), intent(in) ::  vn2,vn3
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(-1:,-1:,-1:), intent(in) ::  B1
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  ns,Ts
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: flagdirich
    real(wp), dimension(:,:,:), intent(in) :: E02src,E03src
    !! these are BG fields use to compute potential source terms
    !! viz. they need to be zeroed out if there is a lagrangian grid...
    real(wp), dimension(:,:), intent(inout) :: Vminx1slab,Vmaxx1slab
    !! need to be able to convert into potential normal deriv.
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: E1,E2,E3,J1,J2,J3
    !! intent(out)
  end subroutine potential_workers_mpi
end interface

interface !< potential_root.f90
  module subroutine potential_root_mpi_curv(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav, &
                                              muP,muH, &
                                              incap,vs2,vs3,vn2,vn3,cfg,B1,ns,Ts,x, &
                                              flagdirich,E02src,E03src,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                                              Vminx1slab,Vmaxx1slab, &
                                              E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)
    integer, intent(in) :: it
    real(wp), intent(in) :: t,dt
    real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
    real(wp), dimension(:,:,:), intent(in) ::  incap
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3
    real(wp), dimension(:,:,:), intent(in) ::  vn2,vn3
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(-1:,-1:,-1:), intent(in) ::  B1
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  ns,Ts
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: flagdirich
    real(wp), dimension(:,:,:), intent(in) :: E02src,E03src
    real(wp), dimension(:,:), intent(inout) :: Vminx1,Vmaxx1    !need to be able to convert these into pot. normal deriv.
    real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    real(wp), dimension(:,:), intent(inout) :: Vminx1slab,Vmaxx1slab    !need to be able to convert into pot. normal deriv.
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: E1,E2,E3,J1,J2,J3
    !! intent(out)
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phiall
    !! not good form, but I'm lazy...  Forgot what I meant by this...
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
  end subroutine potential_root_mpi_curv
end interface

contains
  subroutine electrodynamics_curv(it,t,dt,nn,vn2,vn3,Tn,cfg,ns,Ts,vs1,B1,vs2,vs3,x,efield, &
                           E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec,intvars)
    !! THIS IS A WRAPPER FUNCTION FOR THE ELECTRODYANMICS
    !! PART OF THE MODEL.  BOTH THE ROOT AND WORKER PROCESSES
    !! CALL THIS SAME SUBROUTINE, WHEN THEN BRANCHES INTO
    !! DIFFERENT TASKS FOR EACH AFTER ALL COMPUTE CONDUCTIVITIES
    !! AND INERTIAL CAPACITANCE.
    !!
    !! NOTE THAT THE ALLOCATION STATUS
    !! OF THE ALL VARIABLES FOR THE WORKERS WILL BE UNALLOCATED.
    !! THIS code requires FORTRAN >= 2003 STANDARD.
    integer, intent(in) :: it
    real(wp), intent(in) :: t,dt
    real(wp), dimension(:,:,:,:), intent(in) :: nn
    real(wp), dimension(:,:,:), intent(in) :: vn2,vn3,Tn
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  vs2,vs3
    class(curvmesh), intent(in) :: x
    type(efielddata), intent(inout) :: efield
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: E1,E2,E3,J1,J2,J3
    !! intent(out)
    real(wp), dimension(:,:,:), pointer, intent(inout) :: Phiall
    !! inout since it may not be allocated or deallocated in this procedure
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: sig0,sigPgrav,sigHgrav
    real(wp), dimension(1:lx1,1:lx2,1:lx3,1:lsp) :: muP,muH,nusn
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: incap
    real(wp), dimension(:,:,:), pointer :: sigP,sigH
    real(wp) :: tstart,tfin
    real(wp) :: minh1,maxh1,minh2,maxh2,minh3,maxh3
    ! background variables and boundary conditions, full grid sized variables
    real(wp), dimension(1:x%lx2all,1:x%lx3all), target :: Vminx1,Vmaxx1     !allow pointer aliases for these vars.
    real(wp), dimension(1:x%lx1,1:x%lx3all) :: Vminx2,Vmaxx2
    real(wp), dimension(1:x%lx1,1:x%lx2all) :: Vminx3,Vmaxx3
    ! slab-sized background variables
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: E01,E02,E03,E02src,E03src
    integer :: flagdirich
    real(wp), dimension(1:lx2,1:lx3) :: Vminx1slab,Vmaxx1slab
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: sigNCP,sigNCH


    !> convenience aliases
    sigP=>intvars%sigP; sigH=>intvars%sigH    ! this is the "global" version of conductivity

    !> update conductivities and mobilities
    call cpu_time(tstart)
    call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    if (cfg%flagFBI>1) then     ! need to accumulate nonlinear conductivities if user specifies
      call NLConductivity(nn,Tn,ns,Ts,E2,E3,x,sigP,sigH,sigNCP,sigNCH)
      sigP=sigP+sigNCP
      sigH=sigH+sigNCH
    end if
    call cpu_time(tfin)
    if (mpi_cfg%myid==0) then
      if (cfg%flagcap/=0) then
        if (debug) print *, 'Conductivities and capacitance for time step:  ',t,' took ',tfin-tstart,' seconds...'
      else
        if (debug) print *, 'Conductivities for time step:  ',t,' took ',tfin-tstart,' seconds...'
      end if
    end if

    !> error checking for cap. vs. grid types - viz. we do not support capacitive solves on anything other than Cartesian grids
    if (cfg%flagcap/=0) then      ! some sort of capacitance is being used in the simulation
      call capacitance(ns,B1,cfg,incap)    !> full cfg needed for optional inputs...
      if (it==1) then     !check that we don't have an unsupported grid type for doing capacitance
        minh1=minval(x%h1); maxh1=maxval(x%h1);
        minh2=minval(x%h2); maxh2=maxval(x%h2);
        minh3=minval(x%h3); maxh3=maxval(x%h3);
        if (minh1<0.99d0 .or. maxh1>1.01d0 .or. minh2<0.99d0 .or. maxh2>1.01d0 .or. minh3<0.99d0 .or. maxh3>1.01d0) then
          error stop 'Capacitance is being calculated for possibly unsupported grid type. Please check input file settings.'
        end if
      end if
    else
      incap=0d0
    end if

    !> assign values to the background fields, etc., irrespective of whether or not we do a potential solve
    if (mpi_cfg%myid/=0) then
      call BGfields_boundaries_worker(flagdirich,E01,E02,E03,Vminx1slab,Vmaxx1slab)
    else
      call BGfields_boundaries_root(dt,t,ymd,UTsec,cfg,x,efield, &
                                       flagdirich,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                                       E01,E02,E03,Vminx1slab,Vmaxx1slab)
    end if

    !> must set these variables regardless of whether the solve is done because they are added to the field later
    if (cfg%flaglagrangian) then     ! Lagrangian grid, omit background fields from source terms, note this means that the winds have also been tweaked so that currents/potential source terms will still be correctly computed
      E02src=0._wp; E03src=0._wp
    else                             ! Eulerian grid, use background fields
      E02src=E02; E03src=E03
    end if

    !> Now solve specifically for the *disturbance* potential.  The background electric field will be included in returned electric field and current density values
    if (cfg%potsolve == 1 .or. cfg%potsolve == 3) then    !electrostatic solve or electrostatic alt. solve
      call cpu_time(tstart)
      if (mpi_cfg%myid/=0) then
        !! role-specific communication pattern (all-to-root-to-all), workers initiate with sends
        call potential_workers_mpi(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav,muP,muH,incap,vs2,vs3, &
                                     vn2,vn3,cfg,B1,ns,Ts,x,flagdirich,E02src,E03src, &
                                     Vminx1slab,Vmaxx1slab, &
                                     E1,E2,E3,J1,J2,J3)
      else
        call potential_root_mpi(it,t,dt,sig0,sigP,sigH,sigPgrav,sigHgrav,muP,muH,incap,vs2,vs3,vn2,vn3,cfg,B1,ns,Ts,x, &
                                  flagdirich,E02src,E03src,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                                  Vminx1slab,Vmaxx1slab, &
                                  E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)
      end if
      call cpu_time(tfin)

      if (mpi_cfg%myid==0) then
        if (debug) print *, 'Potential solution for time step:  ',t,' took ',tfin-tstart,' seconds...'
      end if
    else if (cfg%potsolve == 2) then  !inductive form of model, could this be subcycled to speed things up?
      error stop 'Inductive solves are not supported yet; need funding...'
    else   !null solve; set all disturbance fields and currents to zero; these will be accumulated later if backgrounds are used
      E1=0._wp; E2=0._wp; E3=0._wp; J1=0._wp; J2=0._wp; J3=0._wp;
    end if

    !> update *total* electric field variable to include background values
    !if (.not. cfg%flaglagrangian) then     !only add these in if we are not using a lagrangian grid
    !  E2=E2+E02
    !  E3=E3+E03
    !end if

    !> E0?src has been already adjust to account for lagrangian; no need to check again just add
    E2(1:lx1,1:lx2,1:lx3)=E2(1:lx1,1:lx2,1:lx3)+E02src
    E3(1:lx1,1:lx2,1:lx3)=E3(1:lx1,1:lx2,1:lx3)+E03src

    !> velocities should be computed irrespective of whether a solve was done
    call velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
    if (mpi_cfg%myid==0) then
      if (debug) print *, 'Min and max root drift values:  ',minval(vs2),maxval(vs2), minval(vs3),maxval(vs3)
    end if
  end subroutine electrodynamics_curv


  subroutine BGfields_boundaries_root(dt,t,ymd,UTsec,cfg,x,efield, &
                                        flagdirich,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                                        E01,E02,E03,Vminx1slab,Vmaxx1slab)
    !> Root only!  populate arrays for background electric fields and potential/FAC boundary conditions.  If we are root this
    !  involves reading in data from a file or calling a subroutine to fill the arrays and then sending pieces of data
    !  to workers.
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    type(efielddata), intent(inout) :: efield
    integer, intent(out) :: flagdirich
    real(wp), dimension(:,:), intent(inout) :: Vminx1,Vmaxx1
    !! intent(out)
    !! allow pointer aliases for these vars, as required for subroutines internal to potentialBCs.
    real(wp), dimension(:,:), intent(inout) :: Vminx2,Vmaxx2
    !! intent(out)
    real(wp), dimension(:,:), intent(inout)  :: Vminx3,Vmaxx3
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: E01,E02,E03
    !! intent(out)
    real(wp), dimension(:,:), intent(inout) :: Vminx1slab,Vmaxx1slab
    !! intent(out)
    !local work arrays
    real(wp), dimension(:,:,:), allocatable :: E01all,E02all,E03all     !full grid values not needed by root which will collect a source term computation from all workers...
    integer :: iid
    real(wp) :: tstart,tfin
    real(wp), dimension(1:size(Vminx1,1),1:size(Vminx1,2)) :: Vminx1buf,Vmaxx1buf


    !> either read in the data from a file or use a subroutine to set the array values
    allocate(E01all(lx1,lx2all,lx3all),E02all(lx1,lx2all,lx3all),E03all(lx1,lx2all,lx3all))
    call cpu_time(tstart)
    if (cfg%flagE0file==1) then
      call potentialBCs2D_fileinput(dt,t,ymd,UTsec,cfg,x,efield,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                          E01all,E02all,E03all,flagdirich)
    else
      call potentialBCs2D(UTsec,cfg,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                          E01all,E02all,E03all,flagdirich)     !user needs to manually swap x2 and x3 in this function, e.g. for EIA, etc.
    end if
    call cpu_time(tfin)
    if (debug) print *, 'Root has computed BCs in time:  ',tfin-tstart

    do iid=1,mpi_cfg%lid-1
      !! communicate intent for solve to workers so they know whether or not to call mumps fn.
      call mpi_send(flagdirich,1,MPI_INTEGER,iid,tag%flagdirich,MPI_COMM_WORLD)
    end do
    if (debug) print *, 'Root has communicated type of solve to workers:  ',flagdirich

    ! Need to broadcast background fields from root
    ! Need to also broadcast x1 boundary conditions for source term calculations.
    ! This duplicates some code for get_BGEfields, but that is necessary since that lower level routine renormalizes fields using metric factors, which does not need to be done again following a call to potentialBCs
    call bcast_send(E01all,tag%E01,E01)
    call bcast_send(E02all,tag%E02,E02)
    call bcast_send(E03all,tag%E03,E03)
    deallocate(E01all,E02all,E03all)

    !These are pointer targets so don't assume contiguous in memory - pack them into a buffer to be safe
    Vmaxx1buf=Vmaxx1; Vminx1buf=Vminx1;
    call bcast_send(Vminx1buf,tag%Vminx1,Vminx1slab)
    call bcast_send(Vmaxx1buf,tag%Vmaxx1,Vmaxx1slab)
  end subroutine BGfields_boundaries_root


  subroutine BGfields_boundaries_worker(flagdirich,E01,E02,E03,Vminx1slab,Vmaxx1slab)
    !> Worker only!  receive background and boundary condition information from root
    integer, intent(out) :: flagdirich
    real(wp), dimension(:,:,:), intent(inout) :: E01,E02,E03
    !! intent(out)
    real(wp), dimension(:,:), intent(inout) :: Vminx1slab,Vmaxx1slab
    !! intent(out)

    call mpi_recv(flagdirich,1,MPI_INTEGER,0,tag%flagdirich,MPI_COMM_WORLD,MPI_STATUS_IGNORE)

    !Need to broadcast background fields from root
    !Need to also broadcast x1 boundary conditions for source term calculations.
    call bcast_recv(E01,tag%E01)
    call bcast_recv(E02,tag%E02)
    call bcast_recv(E03,tag%E03)
    call bcast_recv(Vminx1slab,tag%Vminx1)
    call bcast_recv(Vmaxx1slab,tag%Vmaxx1)
  end subroutine BGfields_boundaries_worker


  subroutine velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,flaggravdrift,flagdiamagnetic,vs2,vs3)
    !> compute steady state drifts resulting from a range of forces.  Can be used
    !   by both root and worker processes
    real(wp), dimension(:,:,:,:), intent(in) :: muP,muH,nusn
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: E2,E3
    real(wp), dimension(:,:,:), intent(in) :: vn2,vn3
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
    !! these must have ghost cells
    class(curvmesh), intent(in) :: x
    logical, intent(in) :: flaggravdrift
    logical, intent(in) :: flagdiamagnetic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
    !! intent(out)
    !! these have ghost cells
    integer :: isp
    real(wp), dimension(-1:lx1+2,-1:lx2+2,-1:lx3+2) :: pressure    ! temp space for computing these
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: gradlp2,gradlp3

    !! FIXME:  Do we really need separate wind mobility or
    !! can we just compute off electrical mobility as done with gravity.
    !! This is necessary because we are not storing the collision frequencies...
    !! For adding pressure terms we may want to just use collision freq.
    !! to avoid a bunch of slightly different mobility arrays.
    !> electric field and wind terms for ion drifts
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2(1:lx1,1:lx2,1:lx3)-muH(:,:,:,isp)*E3(1:lx1,1:lx2,1:lx3)+ &
                        (muP(:,:,:,isp)*vn2-muH(:,:,:,isp)*vn3)*(ms(isp)*nusn(:,:,:,isp)/qs(isp))
      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2(1:lx1,1:lx2,1:lx3)+muP(:,:,:,isp)*E3(1:lx1,1:lx2,1:lx3)+ &
                        (muH(:,:,:,isp)*vn2+muP(:,:,:,isp)*vn3)*ms(isp)*nusn(:,:,:,isp)/qs(isp)
    end do

    !> Pressure/diamagnetic terms (if required)
    if (flagdiamagnetic) then
      do isp=1,lsp
        !> this behaves better when we take the gradient of log pressure
        pressure(1:lx1,1:lx2,1:lx3)=log(ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp))
        call halo_pot(pressure,tag%pressure,x%flagper,.false.)
        gradlp2=grad3D2(pressure(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
        gradlp3=grad3D3(pressure(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
        vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp) &
                   -muP(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp2(1:lx1,1:lx2,1:lx3) &
                   +muH(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp3(1:lx1,1:lx2,1:lx3)
        vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp) &
                   -muH(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp2(1:lx1,1:lx2,1:lx3) &
                   -muP(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp3(1:lx1,1:lx2,1:lx3)
      end do
    end if

    !> Gravitational drift terms (if required)
    if (flaggravdrift) then
      do isp=1,lsp
        vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(muP(:,:,:,isp)*x%g2-muH(:,:,:,isp)*x%g3)    !FIXME: +muH looks suspicious, I'm changing to (-)
        vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(muH(:,:,:,isp)*x%g2+muP(:,:,:,isp)*x%g3)
      end do
    end if


    !! If it were appropriate this is how polarzations drifts could be computed.  However the particular quasistatic
    !   model that we use explicitly omits this from the drift calculation which is then used in convective term in
    !   polarization current.  Physically it accounts for charge accumulation from polarization currents but not for
    !   the *direct* effect of the polarization term on drift.
    !    do isp=1,lsp
           !! To leading order the ion drifts do not include the polarization parts,
           !! otherwise it may mess up polarization convective term in the electrodynamics solver...
    !      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE2Dt
    !      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE3Dt
    !    end do
  end subroutine velocities


  subroutine potential_sourceterms(sigP,sigH,sigPgrav,sigHgrav,E02,E03,vn2,vn3,B1,muP,muH,ns,Ts,x, &
                                   flaggravdrift,flagdiamagnetic,flagnodivJ0,srcterm)
    !> Compute source terms (inhomogeneous terms) for the potential equation to be solved.  Both root and workers
    !   should be able to use this routine
    real(wp), dimension(:,:,:), intent(in) :: sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:), intent(in) :: E02,E03,vn2,vn3
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
    !! ghost cells
    real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
    !! ghost cells
    class(curvmesh), intent(in) :: x
    logical, intent(in) :: flaggravdrift
    logical, intent(in) :: flagdiamagnetic
    logical, intent(in) :: flagnodivJ0
    real(wp), dimension(:,:,:), intent(inout) :: srcterm
    !! intent(out)
    real(wp), dimension(-1:lx1+2,-1:lx2+2,-1:lx3+2) :: J1,J2,J3    ! why are these local?!
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: divtmp
    !! one extra grid point on either end to facilitate derivatives
    !! haloing assumes existence of two ghost cells

    !-------
    !CONDUCTION CURRENT BACKGROUND SOURCE TERMS FOR POTENTIAL EQUATION. MUST COME AFTER CALL TO BC CODE.
    J1 = 0
    !! so this div is only perp components
    J2 = 0
    J3 = 0
    !! zero everything out to initialize since *accumulating* sources
    if (.not. flagnodivJ0) then
      call acc_perpBGconductioncurrents(sigP,sigH,E02,E03,J2,J3)     !background conduction currents only
      if (debug .and. mpi_cfg%myid==0) print *, 'Workers have computed background field currents...'
    end if
    call acc_perpwindcurrents(sigP,sigH,vn2,vn3,B1,J2,J3)     ! always include wind effects
    if (debug .and. mpi_cfg%myid==0) print *, 'Workers have computed wind currents...'
    if (flagdiamagnetic) then
      call acc_pressurecurrents(muP,muH,ns,Ts,x,J2,J3)
      if (debug .and. mpi_cfg%myid==0) print *, 'Workers have computed pressure currents...'
    end if
    if (flaggravdrift) then
      call acc_perpgravcurrents(sigPgrav,sigHgrav,x%g2,x%g3,J2,J3)
      if (debug .and. mpi_cfg%myid==0) print *, 'Workers have computed gravitational currents...'
    end if

    call halo_pot(J1,tag%J1,x%flagper,.false.)
    call halo_pot(J2,tag%J2,x%flagper,.false.)
    call halo_pot(J3,tag%J3,x%flagper,.false.)

    divtmp=div3D(J1(0:lx1+1,0:lx2+1,0:lx3+1),J2(0:lx1+1,0:lx2+1,0:lx3+1), &
                 J3(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
    srcterm=divtmp(1:lx1,1:lx2,1:lx3)
    !-------
  end subroutine potential_sourceterms


  !> only to be used with electric field arrays that do not have ghost cells
  subroutine acc_perpBGconductioncurrents(sigP,sigH,E2,E3,J2,J3)
    !> ***Accumulate*** conduction currents into the variables J2,J3.  This
    !    routine will not independently add background fields unless they are
    !    already included in E2,3.  The currents are inout meaning that they
    !    must be initialized to zero if you want only the conduction currents,
    !    otherwise this routine just adds to whatever is already in J2,3.
    real(wp), dimension(:,:,:), intent(in) :: sigP,sigH
    real(wp), dimension(1:,1:,1:), intent(in) :: E2,E3
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J2, J3

    J2(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)+sigP*E2(1:lx1,1:lx2,1:lx3)-sigH*E3(1:lx1,1:lx2,1:lx3)
    J3(1:lx1,1:lx2,1:lx3)=J3(1:lx1,1:lx2,1:lx3)+sigH*E2(1:lx1,1:lx2,1:lx3)+sigP*E3(1:lx1,1:lx2,1:lx3)
  end subroutine acc_perpBGconductioncurrents


  subroutine acc_perpconductioncurrents(sigP,sigH,E2,E3,J2,J3)
    !> ***Accumulate*** conduction currents into the variables J2,J3.  This
    !    routine will not independently add background fields unless they are
    !    already included in E2,3.  The currents are inout meaning that they
    !    must be initialized to zero if you want only the conduction currents,
    !    otherwise this routine just adds to whatever is already in J2,3.
    real(wp), dimension(:,:,:), intent(in) :: sigP,sigH
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: E2,E3    ! if used with background could have different lbound so don't assume -1
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J2, J3

    J2(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)+sigP*E2(1:lx1,1:lx2,1:lx3)-sigH*E3(1:lx1,1:lx2,1:lx3)
    J3(1:lx1,1:lx2,1:lx3)=J3(1:lx1,1:lx2,1:lx3)+sigH*E2(1:lx1,1:lx2,1:lx3)+sigP*E3(1:lx1,1:lx2,1:lx3)
  end subroutine acc_perpconductioncurrents


  subroutine acc_perpwindcurrents(sigP,sigH,vn2,vn3,B1,J2,J3)
    !> ***Accumulate*** wind currents into the variables J2,J3.  See conduction currents
    !    routine for additional caveats.
    real(wp), dimension(:,:,:), intent(in) :: sigP,sigH
    real(wp), dimension(:,:,:), intent(in) :: vn2,vn3
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J2, J3

    !! FIXME:  signs here require some explanation...  Perhaps add to formulation doc?
    J2(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)+sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+ &
                            sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
    J3(1:lx1,1:lx2,1:lx3)=J3(1:lx1,1:lx2,1:lx3)+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)- &
                            sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
  end subroutine acc_perpwindcurrents


  subroutine acc_pressurecurrents(muP,muH,ns,Ts,x,J2,J3)
    !> ***Accumulate*** pressure currents into the variables J2,J3.  See conduction currents
    !    routine for additional caveats.
    real(wp), dimension(:,:,:,:), intent(in) :: muP,muH
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
    class(curvmesh), intent(in) :: x
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J2, J3
    integer :: isp
    real(wp), dimension(-1:lx1+2,-1:lx2+2,-1:lx3+2) :: pressure
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: gradp2,gradp3

    do isp=1,lsp
      pressure(1:lx1,1:lx2,1:lx3)=ns(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)
      call halo_pot(pressure,tag%pressure,x%flagper,.false.)
      gradp2=grad3D2(pressure(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
      gradp3=grad3D3(pressure(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
      J2(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)-muP(:,:,:,isp)*gradp2(1:lx1,1:lx2,1:lx3)+ &
                               muH(:,:,:,isp)*gradp3(1:lx1,1:lx2,1:lx3)
      J3(1:lx1,1:lx2,1:lx3)=J3(1:lx1,1:lx2,1:lx3)-muH(:,:,:,isp)*gradp2(1:lx1,1:lx2,1:lx3)- &
                               muP(:,:,:,isp)*gradp3(1:lx1,1:lx2,1:lx3)
    end do
  end subroutine acc_pressurecurrents


  subroutine acc_perpgravcurrents(sigPgrav,sigHgrav,g2,g3,J2,J3)
    !> ***Accumulate*** gravitational currents into the variables J2,J3.  See conduction currents
    !    routine for additional caveats.
    real(wp), dimension(:,:,:), intent(in) :: sigPgrav,sigHgrav
    real(wp), dimension(:,:,:), intent(in) :: g2,g3
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J2, J3

    J2(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)+sigPgrav*g2-sigHgrav*g3
    J3(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)+sigHgrav*g2+sigPgrav*g3
  end subroutine acc_perpgravcurrents


  subroutine pot2perpfield(Phi,x,E2,E3)
    !> computes electric field (perp components only) from a worker potential pattern.  Can
    !   be called by either root or worker processes
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phi
    class(curvmesh), intent(in) :: x
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: E2,E3
    !! intent(out)
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: gradtmp
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: Phitmp
    !! one extra grid point on either end to facilitate derivatives
    !! haloing assumes existence of two ghost cells

    !CALCULATE PERP FIELDS FROM POTENTIAL
    !      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))
    !! causes major memory leak. maybe from arithmetic statement argument?
    !! Left here as a 'lesson learned' (or is it a gfortran bug...)
    !      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))

    call halo_pot(Phi,tag%J1,x%flagper,.true.)
    !call halo_pot(Phihalo,tag%J1,x%flagper,.false.)

    Phitmp=Phi(0:lx1+1,0:lx2+1,0:lx3+1)
    gradtmp=grad3D2(Phitmp,x,0,lx1+1,0,lx2+1,0,lx3+1)   ! FIXME: don't need copy of array???
    E2(1:lx1,1:lx2,1:lx3)=-1*gradtmp(1:lx1,1:lx2,1:lx3)
    gradtmp=grad3D3(Phitmp,x,0,lx1+1,0,lx2+1,0,lx3+1)
    E3(1:lx1,1:lx2,1:lx3)=-1*gradtmp(1:lx1,1:lx2,1:lx3)
    !--------
  end subroutine pot2perpfield


  subroutine parallel_currents(cfg,x,J2,J3,Vminx1slab,Vmaxx1slab,Phi,sig0,flagdirich,J1,E1)
    !> Compute the parallel currents given a potential solution and calculation of perpendicular currents
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J2,J3
    real(wp), dimension(:,:), intent(in) :: Vminx1slab,Vmaxx1slab
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: Phi
    real(wp), dimension(:,:,:), intent(in) :: sig0
    integer, intent(in) :: flagdirich
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J1,E1
    !! intent(out)
    !> local work arrays
    integer :: ix1
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: divtmp
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: divJperp
    !! haloing assumes existence of two ghost cells

    !> inputs to this block of code:  cfg,x,J2,3,Vmaxx1slab,Vminx1slab,flagdirich,gridflag
    !> outputs to code:  J1
    !NOW DEAL WITH THE PARALLEL FIELDS AND ALL CURRENTS
    if (lx2/=1 .and. lx3/=1 .and. cfg%potsolve ==1) then    !we did a field-integrated solve for potential
      if (debug) print*, 'Appear to need to differentiate to get J1...'

      !-------
      !NOTE THAT A DIRECT E1ALL CALCULATION WILL GIVE ZERO, SO USE INDIRECT METHOD, AS FOLLOWS
      J1= 0    !a placeholder so that only the perp divergence is calculated - will get overwritten later.
    !      divJperp=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)

      if (cfg%flagJpar) then   ! user can elect not to compute Jpar, which can be prone to artifacts particularly at low resolution
        call halo_pot(J1,tag%J1,x%flagper,.false.)
        call halo_pot(J2,tag%J2,x%flagper,.false.)
        call halo_pot(J3,tag%J3,x%flagper,.false.)

        divtmp=div3D(J1(0:lx1+1,0:lx2+1,0:lx3+1),J2(0:lx1+1,0:lx2+1,0:lx3+1), &
                     J3(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
        divJperp=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*divtmp(1:lx1,1:lx2,1:lx3)
        if (flagdirich /= 1) then
          !! Neumann conditions, this is boundary location-agnostic since both bottom and top FACs are known
          !! - they have to  be loaded into VVmaxx1 and Vminx1
    !      if (debug) print *, 'Nuemann boundaries, integrating from highest altitude down to preserve accuracy...'
          if (gridflag==0) then     !closed dipole grid, really would be best off integrating from the source hemisphere
    !        if (debug) print *,  'Closed dipole grid; integration starting in source hemisphere (if applicable)...', &
    !                       minval(Vmaxx1slab), &
    !                       maxval(Vmaxx1slab)
            if (cfg%sourcemlat >= 0) then    !integrate from northern hemisphere
    !          if (debug) print *, 'Source is in northern hemisphere (or there is no source)...'
              J1(1:lx1,1:lx2,1:lx3)=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
              do ix1=1,lx1
                J1(ix1,1:lx2,1:lx3)= 1/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                                 (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,1:lx2,1:lx3))
              end do
            else
    !          if (debug) print *, 'Source in southern hemisphere...'
              J1(1:lx1,1:lx2,1:lx3)=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current starting from minx1
              do ix1=1,lx1
                J1(ix1,1:lx2,1:lx3)= 1/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                                 (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,1:lx2,1:lx3))
              end do
            end if
          ! For open grids always integrate from the bottom (zero current)
          elseif (gridflag==1) then    !this would be an inverted grid, this max altitude corresponds to the min value of x1
    !        if (debug) print *,  'Inverted grid; integration starting at min x1 (highest alt. or southern hemisphere)...', &
    !                       minval(Vminx1slab), &
    !                       maxval(Vminx1slab)
            J1(1:lx1,1:lx2,1:lx3)=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
            do ix1=1,lx1
              J1(ix1,1:lx2,1:lx3)= 1/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                               (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,1:lx2,1:lx3))
            end do

          else        !minx1 is at the bottom of the grid to integrate from max x1
    !        if (debug) print *,  'Non-inverted grid; integration starting at max x1...', minval(Vmaxx1slab), maxval(Vmaxx1slab)
            J1(1:lx1,1:lx2,1:lx3)=integral3D1_curv_alt(divJperp,x,1,lx1)
            !! int divperp of BG current, go from maxval(x1) to location of interest
            do ix1=1,lx1
              J1(ix1,1:lx2,1:lx3)= 1/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                               (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,1:lx2,1:lx3))
            end do

          end if
        else
          !! Dirichlet conditions - we need to integrate from the ***lowest altitude***
          !! (where FAC is known to be zero, note this is not necessarilty the logical bottom of the grid), upwards (to where it isn't)
          if (gridflag/=2) then    !inverted grid (logical top is the lowest altitude)
            if (debug) print *, 'Inverted grid detected - integrating logical top downward to compute FAC...'
            J1(1:lx1,1:lx2,1:lx3)=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current
            do ix1=1,lx1
              J1(ix1,1:lx2,1:lx3)= 1/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                               (J1(ix1,1:lx2,1:lx3))    !FAC AT TOP ASSUMED TO BE ZERO
            end do
          else      !non-inverted grid (logical bottom is the lowest altitude - so integrate normy)
            if (debug) print *, 'Non-inverted grid detected - integrating logical bottom to top to compute FAC...'
            J1(1:lx1,1:lx2,1:lx3)=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
            do ix1=1,lx1
              J1(ix1,1:lx2,1:lx3)= 1/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                               (-J1(ix1,1:lx2,1:lx3))    !FAC AT THE BOTTOM ASSUMED TO BE ZERO
            end do
          end if
        end if
        E1(1:lx1,1:lx2,1:lx3)=J1(1:lx1,1:lx2,1:lx3)/sig0
        !-------
      end if ! flagJpar
    else   !we resolved the field line (either 2D solve or full 3D) so just differentiate normally
      !-------
      E1(1:lx1,1:lx2,1:lx3)=grad3D1(Phi(1:lx1,1:lx2,1:lx3),x,1,lx1,1,lx2,1,lx3)    !no haloing required since x1-derivative
      E1(1:lx1,1:lx2,1:lx3)=-E1(1:lx1,1:lx2,1:lx3)
      J1(1:lx1,1:lx2,1:lx3)=sig0*E1(1:lx1,1:lx2,1:lx3)
      !print*, 'parallel fields:  ',maxval(abs(E1))
      !-------
    end if
  end subroutine parallel_currents


  subroutine polarization_currents(cfg,x,dt,incap,E2,E3,E2prev,E3prev,v2,v3,J1pol,J2pol,J3pol)
    !> Computes the polarization currents resulting from time-dependence and shearing of the plasma
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: E2,E3
    real(wp), dimension(:,:,:), intent(in) :: incap,E2prev,E3prev,v2,v3
    real(wp), dimension(1:lx1,1:lx2,1:lx3), intent(inout) :: J1pol,J2pol,J3pol
    !! intent(out)
    ! internal work arrays
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: DE2Dt,DE3Dt,grad2E,grad3E
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: divtmp

    ! check whether electrodynamics is being used or not
    if (cfg%flagcap/=0) then
      !differentiate E2 in x2
      call halo_pot(E2,tag%J1,x%flagper,.false.)
      divtmp=grad3D2(E2(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
      grad2E=divtmp(1:lx1,1:lx2,1:lx3)
      !differentiate E2 in x3
      divtmp=grad3D3(E2(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
      grad3E=divtmp(1:lx1,1:lx2,1:lx3)
      !compute total derivative in x2
      DE2Dt=(E2(1:lx1,1:lx2,1:lx3)-E2prev)/dt+v2*grad2E+v3*grad3E

      !differentiate E3 in x2
      call halo_pot(E3,tag%J1,x%flagper,.false.)
      divtmp=grad3D2(E3(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
      grad2E=divtmp(1:lx1,1:lx2,1:lx3)
      !differentiate E3 in x3
      divtmp=grad3D3(E3(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
      grad3E=divtmp(1:lx1,1:lx2,1:lx3)
      !x3 total derivative
      DE3Dt=(E3(1:lx1,1:lx2,1:lx3)-E3prev)/dt+v2*grad2E+v3*grad3E

      !convert derivative into polarization current density
      J1pol= 0
      J2pol=incap*DE2Dt
      J3pol=incap*DE3Dt
    else       !pure electrostatic solve was done, zero out everything
      DE2Dt= 0
      DE3Dt= 0
      J1pol= 0
      J2pol= 0
      J3pol= 0
    end if
  end subroutine polarization_currents


  subroutine get_BGEfields(x,E01,E02,E03,efield)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:), intent(inout) :: E01,E02,E03
    !! intent(out)
    type(efielddata), intent(inout) :: efield
    !> routine to pull the background electric fields from the BCs modules and distribute
    !   them to worker processes; called by both root and worker processes
    real(wp), dimension(:,:,:), allocatable :: E01all,E02all,E03all    !space to pull the full dataset out of the module


    if (mpi_cfg%myid==0) then
      allocate(E01all(lx1,lx2all,lx3all),E02all(lx1,lx2all,lx3all),E03all(lx1,lx2all,lx3all))
      E01all=0._wp
      call compute_rootBGEfields(x,E02all,E03all,efield)

      call bcast_send(E01all,tag%E01,E01)
      call bcast_send(E02all,tag%E02,E02)
      call bcast_send(E03all,tag%E03,E03)
      deallocate(E01all,E02all,E03all)
    else
      call bcast_recv(E01,tag%E01)
      call bcast_recv(E02,tag%E02)
      call bcast_recv(E03,tag%E03)
    end if
  end subroutine get_BGEfields


  subroutine halo_pot(parmhalo,tagcurrent,flagper,flagdegrade)
    !THIS SUBROUTINE REPLICATES A COMMON MESSAGE PASSING SCHEME USED IN THE COMPUTATION
    !OF ELECTRODYNAMICS PARAMETERS THAT RESULT FROM DERIVATIVE (WHICH REQUIRE HALOING)
    real(wp), intent(inout), dimension(-1:,-1:,-1:) :: parmhalo
    integer, intent(in) :: tagcurrent
    logical, intent(in) :: flagper
    logical, intent(in) :: flagdegrade    !whether or not to degrade edge derivatives to first order in x2
    integer :: lx1,lx2,lx3
    integer :: idleft,idright,iddown,idup


    lx1=size(parmhalo,1)-4
    lx2=size(parmhalo,2)-4
    lx3=size(parmhalo,3)-4

    idleft=mpi_cfg%myid3-1
    idright=mpi_cfg%myid3+1
    iddown=mpi_cfg%myid2-1
    idup=mpi_cfg%myid2+1

    parmhalo(0,1:lx2,1:lx3)=parmhalo(1,1:lx2,1:lx3)
    parmhalo(lx1+1,1:lx2,1:lx3)=parmhalo(lx1,1:lx2,1:lx3)

    call halo(parmhalo,1,tagcurrent,flagper)     !this particular type of message passing only needs a single ghost cell

    ! x2 global boundary
    if (iddown==-1) then
      if (flagdegrade .and. lx2>1) then     !for whatever reason this fails ctest without checking lx2>1
        parmhalo(1:lx1,0,1:lx3)=-1*parmhalo(1:lx1,2,1:lx3)+2*parmhalo(1:lx1,1,1:lx3)
      else
        parmhalo(1:lx1,0,1:lx3)=parmhalo(1:lx1,1,1:lx3)
      end if
    end if
    if (idup==mpi_cfg%lid2) then
      if (flagdegrade .and. lx2>1) then
        parmhalo(1:lx1,lx2+1,1:lx3)=2*parmhalo(1:lx1,lx2,1:lx3)-parmhalo(1:lx1,lx2-1,1:lx3)
      else
        parmhalo(1:lx1,lx2+1,1:lx3)=parmhalo(1:lx1,lx2,1:lx3)
      end if
    end if

    ! x3 global boundary
    if (.not. flagper) then
      !! mustn't overwrite ghost cells if perioidc is chosen
      if (idleft==-1) then
        parmhalo(1:lx1,1:lx2,0)=parmhalo(1:lx1,1:lx2,1)
      end if
      if (idright==mpi_cfg%lid3) then
        parmhalo(1:lx1,1:lx2,lx3+1)=parmhalo(1:lx1,1:lx2,lx3)
      end if
    end if
  end subroutine halo_pot
end module potential_comm
