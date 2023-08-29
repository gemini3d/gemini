!> mpi-related gemini functionality
module gemini3d_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: iso_c_binding, only : c_ptr
use phys_consts, only: wp,debug
use mpimod, only: mpi_manualgrid, process_grid_auto, mpi_cfg, mpibreakdown, mpisetup, tag=>gemini_mpi, halo
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use io, only: output_plasma,output_aur,find_milestone,input_plasma,create_outdir
use potential_comm, only: get_BGEfields,velocities
use grid, only: lx1,lx2,lx3, grid_drift, read_grid, calc_subgrid_size
use collisions, only: conductivities
use potentialBCs_mumps, only: init_Efieldinput
use potential_comm,only : pot2perpfield, electrodynamics
use neutral_perturbations, only: init_neutralperturb,neutral_denstemp_update,neutral_wind_update,neutral_perturb
use temporal_mpi, only : dt_comm
use sanity_check, only : check_finite_pertub, check_finite_output
use advec_mpi, only: halo_interface_vels_allspec
use multifluid_mpi, only: halo_allparams, halo_fluidvars
use sources_mpi, only: RK2_prep_mpi_allspec, RK2_global_boundary_allspec
use ionization_mpi, only: get_gavg_Tinf
use neutral_perturbations, only: clear_dneu
use gemini3d, only: fluidvar_pointers,fluidauxvar_pointers, electrovar_pointers, gemini_work,  &
                      v2grid, v3grid, setv2v3, set_start_timefromcfg
use gemini_work_def, only: gemini_work

implicit none (type, external)
private
public :: init_procgrid, outdir_fullgridvaralloc, read_grid_in, get_initial_state, &
            BGfield_Lagrangian, check_dryrun, check_fileoutput,  &
            get_initial_drifts, init_Efieldinput_in, pot2perpfield_in, init_neutralperturb_in, dt_select, &
            neutral_atmos_wind_update, neutral_perturb_in, electrodynamics_in, check_finite_output_in, &
            halo_interface_vels_allspec_in, halo_allparams_in, &
            RK2_prep_mpi_allspec_in,get_gavg_Tinf_in, clear_dneu_in, mpisetup_in, mpiparms, calc_subgrid_size_in, &
            RK2_global_boundary_allspec_in, halo_fluidvars_in

real(wp), parameter :: dtscale=2                     ! controls how rapidly the time step is allowed to change

contains
  !> call the mpi setup from our module
  subroutine mpisetup_in()
    call mpisetup()
  end subroutine mpisetup_in


  !> my id and the total number of processes for this run
  subroutine mpiparms(myid,lid)
    integer, intent(inout) :: myid,lid

    myid=mpi_cfg%myid
    lid=mpi_cfg%lid
  end subroutine mpiparms


  !> create output directory and allocate full grid potential storage
  subroutine outdir_fullgridvaralloc(cfg,intvars,lx1,lx2all,lx3all)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), intent(inout) :: intvars
    integer, intent(in) :: lx1,lx2all,lx3all

    !> create a place, if necessary, for output datafiles
    if (mpi_cfg%myid==0) then
      call create_outdir(cfg)
      allocate(intvars%Phiall(-1:lx1+2,-1:lx2all+2,-1:lx3all+2))
    end if
  end subroutine outdir_fullgridvaralloc


  !> read in the grid and distribute to workers
  subroutine read_grid_in(cfg,x)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), pointer, intent(inout) :: x

    call read_grid(cfg%indatsize,cfg%indatgrid,cfg%flagperiodic, x)
    !! read in a previously generated grid from filenames listed in input file
  end subroutine read_grid_in


  !> interface for setting simulation subgrid sizes for a particular worker in grid module
  subroutine calc_subgrid_size_in(lx2all,lx3all)
    integer, intent(in) :: lx2all, lx3all

    call calc_subgrid_size(lx2all,lx3all)
  end subroutine calc_subgrid_size_in


!  !> load initial conditions and check if this is a restart run; set time variables accordingly
!  subroutine get_initial_state(cfg,fluidvars,electrovars,intvars,x,UTsec,ymd,tdur)
!    type(gemini_cfg), intent(inout) :: cfg
!    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
!    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars
!    type(gemini_work), intent(inout) :: intvars
!    class(curvmesh), intent(in) :: x
!    real(wp), intent(inout) :: UTsec
!    integer, dimension(3), intent(inout) :: ymd
!    real(wp), intent(inout) :: tdur
!
!    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
!    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    integer, dimension(3) :: ymdtmp
!    real(wp) :: UTsectmp,ttmp
!    character(:), allocatable :: filetmp
!
!    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!
!    call find_milestone(cfg, ttmp, ymdtmp, UTsectmp, filetmp)
!    if ( ttmp > 0 ) then
!      !! restart scenario
!      if (mpi_cfg%myid==0) then
!        print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!        print*, '! Restarting simulation from time:  ',ymdtmp,UTsectmp
!        print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!      end if
!
!      !! Set start variables accordingly and read in the milestone
!      UTsec=UTsectmp
!      ymd=ymdtmp
!
!      ! FIXME: instead keep tdur and just adjust the start time of the simulation to be closer to endtime
!      tdur=cfg%tdur-ttmp    ! subtract off time that has elapsed to milestone
!      ! FIXME: need to feed in the t variable and overwrite it if we are restarting.  
!
!      if (mpi_cfg%myid==0) then
!        print*, 'Treating the following file as initial conditions:  ',filetmp
!        print*, ' full duration:  ',cfg%tdur,'; remaining simulation time:  ',tdur    ! FIXME: should be tdur-t once t adjusted
!      end if
!
!      if (tdur <= 1e-6_wp .and. mpi_cfg%myid==0) error stop 'Cannot restart simulation from the final time step!'
!
!      cfg%tdur=tdur         ! just to insure consistency
!      call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,filetmp,ns,vs1,Ts,Phi,intvars%Phiall)
!    else !! start at the beginning
!      ! UTsec = cfg%UTsec0
!      ! ymd = cfg%ymd0
!      ! tdur = cfg%tdur
!      call set_start_timefromcfg(cfg,ymd,UTsec,tdur)
!
!      if (tdur <= 1e-6_wp .and. mpi_cfg%myid==0) error stop 'Simulation is of zero time duration'
!      print*, 'Starting from beginning of simulation...'
!      call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,cfg%indatfile,ns,vs1,Ts,Phi,intvars%Phiall)
!    end if
!  end subroutine get_initial_state


  !> load initial conditions and check if this is a restart run; set time variables accordingly
  subroutine get_initial_state(cfg,fluidvars,electrovars,intvars,x,UTsec,ymd,tdur,t)
    type(gemini_cfg), intent(inout) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(inout) :: UTsec
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: tdur,t

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    integer, dimension(3) :: ymdtmp
    real(wp) :: UTsectmp
    character(:), allocatable :: filetmp
    real(wp) :: tremaining

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    call find_milestone(cfg, t, ymdtmp, UTsectmp, filetmp)
    if ( t > 0 ) then
      !! Set start variables accordingly and read in the milestone
      UTsec=UTsectmp
      ymd=ymdtmp

      tremaining=cfg%tdur-t    ! subtract off time that has elapsed to milestone
      tdur=cfg%tdur

      if (mpi_cfg%myid==0) then
        print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print*, '! Restarting simulation from time:  ',ymdtmp,UTsectmp
        print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print*, 'Treating the following file as initial conditions:  ',filetmp
        print*, ' full duration:  ',cfg%tdur,'; remaining simulation time:  ',tremaining
        print*, ' simulation start time:  ',t
        if (tremaining<=1e-6_wp) error stop 'Cannot restart simulation from the final time step!'
      end if

      !cfg%tdur=tdur         ! just to insure consistency
      call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,filetmp,ns,vs1,Ts,Phi,intvars%Phiall)
    else !! start at the beginning
      ! UTsec = cfg%UTsec0
      ! ymd = cfg%ymd0
      ! tdur = cfg%tdur
      call set_start_timefromcfg(cfg,ymd,UTsec,tdur)
      t=0._wp
      tdur=cfg%tdur

      if (tdur <= 1e-6_wp .and. mpi_cfg%myid==0) error stop 'Simulation is of zero time duration'
      print*, 'Starting from beginning of simulation...'
      call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,cfg%indatfile,ns,vs1,Ts,Phi,intvars%Phiall)
    end if
  end subroutine get_initial_state


  !> check whether file output should be done and complete it
  subroutine check_fileoutput(cfg,fluidvars,electrovars,intvars,t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: tout,tglowout,tmilestone
    integer, intent(inout) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp) :: tstart,tfin
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    if (abs(t-tout) < 1d-5) then
      tout = tout + cfg%dtout
      if (cfg%nooutput ) then
        if (mpi_cfg%myid==0) write(stderr,*) 'WARNING: skipping file output at sim time (sec)',t
        return
      endif
      !! close enough to warrant an output now...
      if (mpi_cfg%myid==0 .and. debug) call cpu_time(tstart)

      !! We may need to adjust flagoutput if we are hitting a milestone
      flagoutput=cfg%flagoutput
      if (cfg%mcadence>0 .and. abs(t-tmilestone) < 1d-5) then
        flagoutput=1    !force a full output at the milestone
        call output_plasma(cfg%outdir,flagoutput,ymd, &
          UTsec,vs2,vs3,ns,vs1,Ts,intvars%Phiall,J1,J2,J3, &
          cfg%out_format,intvars)
        tmilestone = t + cfg%dtout * cfg%mcadence
        if(mpi_cfg%myid==0) print*, 'Milestone output triggered.'
      else
        call output_plasma(cfg%outdir,flagoutput,ymd, &
          UTsec,vs2,vs3,ns,vs1,Ts,intvars%Phiall,J1,J2,J3, &
          cfg%out_format,intvars)
      end if
      if (mpi_cfg%myid==0 .and. debug) then
        call cpu_time(tfin)
        print *, 'Plasma output done for time step:  ',t,' in cpu_time of:  ',tfin-tstart
      endif
    end if

    !> GLOW file output
    if ((cfg%flagglow /= 0) .and. (abs(t-tglowout) < 1d-5)) then !same as plasma output
      call cpu_time(tstart)
      call output_aur(cfg%outdir, cfg%flagglow, ymd, UTsec, intvars%iver, cfg%out_format)
      if (mpi_cfg%myid==0) then
        call cpu_time(tfin)
        print *, 'Auroral output done for time step:  ',t,' in cpu_time of: ',tfin-tstart
      end if
      tglowout = tglowout + cfg%dtglowout
    end if
  end subroutine


  !> check whether a dryrun simulation was done
  subroutine check_dryrun(cfg)
    type(gemini_cfg), intent(in) :: cfg
    character(8) :: date
    character(10) :: time
    integer :: ierr

    if (cfg%dryrun) then
      ierr = mpibreakdown()
      if (ierr /= 0) error stop 'Gemini dry run MPI shutdown failure'
      call date_and_time(date,time)
      print '(/,A)', 'DONE: ' // date(1:4) // '-' // date(5:6) // '-' // date(7:8) // 'T' &
        // time(1:2) // ':' // time(3:4) // ':' // time(5:)
      stop "OK: Gemini dry run"
    endif
  end subroutine check_dryrun


  !> prep simulation for use of Lagrangian grid, if needed
  subroutine BGfield_Lagrangian(cfg,x,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    real(wp) :: v2gridtmp,v3gridtmp

    real(wp), dimension(:,:,:), allocatable :: E01,E02,E03
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)     ! bind pointers for EM variables

    allocate(E01(lx1,lx2,lx3),E02(lx1,lx2,lx3),E03(lx1,lx2,lx3))    ! allocate without ghost cells
    E01=0; E02=0; E03=0;
    if (cfg%flagE0file==1) then
      call get_BGEfields(x,E01,E02,E03,intvars%efield)
    end if
    if (cfg%flaglagrangian) then    ! Lagrangian (moving) grid; compute from input background electric fields
      call grid_drift(x,E02,E03,v2gridtmp,v3gridtmp)
      call setv2v3(v2gridtmp,v3gridtmp)
      if (mpi_cfg%myid==0) print*, mpi_cfg%myid,' using Lagrangian grid moving at:  ',v2grid,v3grid
    else                            ! stationary grid
      v2gridtmp=0._wp
      v3gridtmp=0._wp
      call setv2v3(v2gridtmp,v3gridtmp)
      E1(1:lx1,1:lx2,1:lx3) = E1(1:lx1,1:lx2,1:lx3) + E01    ! FIXME: this is before dist fields are computed???
      E2(1:lx1,1:lx2,1:lx3) = E2(1:lx1,1:lx2,1:lx3) + E02
      E3(1:lx1,1:lx2,1:lx3) = E3(1:lx1,1:lx2,1:lx3) + E03
    end if
    deallocate(E01,E02,E03)

    if (mpi_cfg%myid==0) then
      print*, 'Recomputed initial fields including background (BGfieldLagrangian):'
      print*, '    ',minval(E1(1:lx1,1:lx2,1:lx3)),maxval(E1(1:lx1,1:lx2,1:lx3))
      print*, '    ',minval(E2(1:lx1,1:lx2,1:lx3)),maxval(E2(1:lx1,1:lx2,1:lx3))
      print*, '    ',minval(E3(1:lx1,1:lx2,1:lx3)),maxval(E3(1:lx1,1:lx2,1:lx3))
    end if
  end subroutine BGfield_Lagrangian


  !> initial drifts at the start of the simulation
  subroutine get_initial_drifts(cfg,x,fluidvars,fluidauxvars,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(in) :: intvars

    real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:,:), allocatable :: muP,muH,nusn
    integer :: lx1,lx2,lx3,lsp
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    ! calculate drifts
    lx1=x%lx1; lx2=x%lx2; lx3=x%lx3; lsp=size(ns,4);
    allocate(sig0(lx1,lx2,lx3),sigP(lx1,lx2,lx3),sigH(lx1,lx2,lx3),sigPgrav(lx1,lx2,lx3),sigHgrav(lx1,lx2,lx3))
    allocate(muP(lx1,lx2,lx3,lsp),muH(lx1,lx2,lx3,lsp),nusn(lx1,lx2,lx3,lsp))
    call conductivities(intvars%atmos%nn,intvars%atmos%Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    call velocities(muP,muH,nusn,E2,E3,intvars%atmos%vn2,intvars%atmos%vn3,ns,Ts,x, &
                      cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
    deallocate(sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    if(mpi_cfg%myid==0) then
      print*, 'Recomputed initial drifts:  '
      print*, '    ',minval(vs2(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs2(1:lx1,1:lx2,1:lx3,1:lsp))
      print*, '    ',minval(vs3(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs3(1:lx1,1:lx2,1:lx3,1:lsp))
    end if
  end subroutine get_initial_drifts


  !> initialize the process gridf for this simulation
  subroutine init_procgrid(lx2all,lx3all,lid2in,lid3in)
    integer, intent(in) :: lx2all,lx3all,lid2in,lid3in

    if (lid2in==-1) then
      call process_grid_auto(lx2all, lx3all)
      !! grid_size defines lx2all and lx3all
    else
      call mpi_manualgrid(lx2all, lx3all, lid2in, lid3in)
    endif
    print '(A, I0, A1, I0)', 'process grid (Number MPI processes) x2, x3:  ',mpi_cfg%lid2, ' ', mpi_cfg%lid3
    print '(A, I0, A, I0, A1, I0)', 'Process:',mpi_cfg%myid,' at process grid location: ',mpi_cfg%myid2,' ',mpi_cfg%myid3
  end subroutine init_procgrid


  !> initialize electric field input data
  subroutine init_Efieldinput_in(cfg,x,dt,intvars,ymd,UTsec)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    type(gemini_work), intent(inout) :: intvars
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call init_Efieldinput(dt,cfg,ymd,UTsec,x,intvars%efield)
  end subroutine init_Efieldinput_in


  !> convert potential to electric field by differentiating
  subroutine pot2perpfield_in(x,electrovars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars

    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    E1 = 0
    call pot2perpfield(Phi,x,E2,E3)
    if(mpi_cfg%myid==0) then
      print '(A)', 'Recomputed initial dist. fields via pot2perpfield:'
      print*, '    gemini ',minval(E1(1:lx1,1:lx2,1:lx3)),maxval(E1(1:lx1,1:lx2,1:lx3))
      print*, '    gemini ',minval(E2(1:lx1,1:lx2,1:lx3)),maxval(E2(1:lx1,1:lx2,1:lx3))
      print*, '    gemini ',minval(E3(1:lx1,1:lx2,1:lx3)),maxval(E3(1:lx1,1:lx2,1:lx3))
      print*, '    gemini ',minval(Phi(1:lx1,1:lx2,1:lx3)),maxval(Phi(1:lx1,1:lx2,1:lx3))
    end if
  end subroutine pot2perpfield_in


  !> initialize neutral perturbations
  subroutine init_neutralperturb_in(dt,cfg,x,intvars,ymd,UTsec)
    real(wp), intent(in) :: dt
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    type(gemini_work), intent(inout) :: intvars
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call init_neutralperturb(cfg,x,dt,ymd,UTsec,intvars%atmosperturb)
  end subroutine init_neutralperturb_in


  !> select time step and throttle if changing too rapidly
  subroutine dt_select(cfg,x,fluidvars,fluidauxvars,it,t,tout,tglowout,dt)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    integer, intent(in) :: it
    real(wp), intent(in) :: t,tout,tglowout
    real(wp), intent(inout) :: dt

    real(wp) :: dtprev
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    ! bind pointers for fluid and auxiliary variables
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !> save prior time step
    dtprev = dt

    !> time step calculation, requires workers to report their most stringent local stability constraint
    call dt_comm(t,tout,tglowout,cfg,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,dt)

    !> do not allow the time step to change too rapidly
    if (it>1) then
      if(dt/dtprev > dtscale) then
        !! throttle how quickly we allow dt to increase
        dt=dtscale*dtprev
        if (mpi_cfg%myid == 0) then
          print '(A,EN14.3)', 'Throttling dt to:  ',dt
        end if
      end if
    end if
  end subroutine dt_select


  !> apply neutral perturbations/background and assign to main code variables
  subroutine neutral_atmos_wind_update(intvars)
    type(gemini_work), intent(inout) :: intvars

    call neutral_denstemp_update(intvars%atmos,intvars%atmosperturb)
    call neutral_wind_update(v2grid,v3grid,intvars%atmos,intvars%atmosperturb)
  end subroutine neutral_atmos_wind_update


  !> compute neutral perturbations and apply to main code variables
  subroutine neutral_perturb_in(cfg,intvars,x,dt,t,ymd,UTsec)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(inout) :: x     ! unit vectors could be deallocated in this procedure
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    call neutral_perturb(cfg,dt,t,ymd,UTsec,x,v2grid,v3grid,intvars%atmos,intvars%atmosperturb)
    call check_finite_pertub(cfg%outdir, t, mpi_cfg%myid, intvars%atmos%nn, intvars%atmos%Tn, &
                               intvars%atmos%vn1, intvars%atmos%vn2, intvars%atmos%vn3)
  end subroutine neutral_perturb_in


  !> call electrodynamics solution
  subroutine electrodynamics_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,it,t,dt,ymd,UTsec)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: it
    real(wp), intent(in) :: t,dt
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    ! E&M solves
    call electrodynamics(it,t,dt,intvars%atmos%nn,intvars%atmos%vn2,intvars%atmos%vn3,intvars%atmos%Tn, &
                           cfg,ns,Ts,vs1,B1,vs2,vs3,x,intvars%efield,E1,E2,E3,J1,J2,J3, &
                           intvars%Phiall,ymd,UTsec,intvars)
  end subroutine electrodynamics_in


  !> check main state variables for finiteness
  subroutine check_finite_output_in(cfg,fluidvars,electrovars,t)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    real(wp), intent(in) :: t

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    call check_finite_output(cfg%outdir,t,mpi_cfg%myid,vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)
  end subroutine check_finite_output_in

  ! FIXME:  deprecated; is easier/better just to halo all params prior to advection and interface vels
  !           calculation...
  !> haloing for computing cell interface velocities
  subroutine halo_interface_vels_allspec_in(x,fluidvars,lsp)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    integer, intent(in) :: lsp

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call halo_interface_vels_allspec(x%flagper,vs2,vs3,lsp)
  end subroutine halo_interface_vels_allspec_in


  !> halo all ***advected*** parameters
  subroutine halo_allparams_in(x,fluidvars,fluidauxvars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    ! halo all fluid parameters
    call halo_allparams(ns,rhovs1,rhoes,x%flagper)
  end subroutine halo_allparams_in


  !> halo all parameters, including velocities
   subroutine halo_fluidvars_in(x,fluidvars,fluidauxvars)
     class(curvmesh), intent(in) :: x
     real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
     real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars

     real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
     real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
     real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

     ! bind pointers
     call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
     call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

     ! halo all fluid parameters
     call halo_fluidvars(ns,rhovs1,rhoes,vs2,vs3,x%flagper)
   end subroutine halo_fluidvars_in


  !> prepare and then halo data for compression substep
   subroutine RK2_prep_mpi_allspec_in(x,fluidvars)
     class(curvmesh), intent(in) :: x
     real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
     real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
     real(wp), dimension(-1:x%lx1+2,-1:x%lx2+2,-1:x%lx3+2) :: param
     integer :: isp,lsp

     call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
     lsp=size(ns,4)
     call RK2_prep_mpi_allspec(vs1,vs2,vs3,x%flagper)
   end subroutine RK2_prep_mpi_allspec_in


  !> prepare data for compression substep deal with global boundaries
   subroutine RK2_global_boundary_allspec_in(x,fluidvars)
     class(curvmesh), intent(in) :: x
     real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars

     real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

     call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
     call RK2_global_boundary_allspec(vs1,vs2,vs3,x%flagper)
   end subroutine RK2_global_boundary_allspec_in


  !> agree on average value of gravity and exospheric temp
  subroutine get_gavg_Tinf_in(intvars,gavg,Tninf)
    type(gemini_work), intent(in) :: intvars
    real(wp), intent(inout) :: gavg,Tninf

    call get_gavg_Tinf(intvars%atmos,gavg,Tninf)
  end subroutine get_gavg_Tinf_in


  !> deallocate module storage for neutral perturbations
  subroutine clear_dneu_in(intvars)
    type(gemini_work), intent(inout) :: intvars

    call clear_dneu(intvars%atmosperturb)
  end subroutine clear_dneu_in
end module gemini3d_mpi
