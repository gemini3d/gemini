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

!> This module is intended to have various interfaces/wrappers for main gemini functionality
!!   that does not involve mpi or mpi-dependent modules.  There will be a separate module with the
!!   C bindings and pointer conversions that can be called from C (i.e. wrappers for these routines).
!!   For the most part this is a bunch of "getter" routines.
module gemini3d

use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float, c_loc, c_null_ptr, c_ptr, c_f_pointer
use gemini_cli, only : cli
use gemini_init, only : find_config, check_input_files
use phys_consts, only: wp,debug,lnchem,lwave,lsp,pi
use meshobj, only: curvmesh
use precipdataobj, only: precipdata
use efielddataobj, only: efielddata
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use neutraldata3Dobj_fclaw, only: neutraldata3D_fclaw
use gemini3d_config, only: gemini_cfg
use collisions, only: conductivities
use filesystem, only : expanduser
use temporal, only: cflcalc
use grid, only: grid_size,lx1,lx2,lx3,lx2all,lx3all,grid_from_extents,read_size_gridcenter, get_gridcenter, &
                  grid_internaldata_ungenerate, meshobj_alloc, meshobj_dealloc, grid_internaldata_alloc, &
                  grid_internaldata_generate, get_x1coords, get_fullgrid_lims
use gemini3d_config, only : gemini_cfg,read_configfile
use precipBCs_mod, only: init_precipinput
use msis_interface, only : msisinit
use neutral, only: neutral_info,init_neutralBG,neutral_atmos,neutral_winds,neutral_info_alloc,neutral_info_dealloc
use multifluid, only : sweep3_allparams,sweep1_allparams,sweep2_allparams,source_loss_allparams,VNRicht_artvisc,compression, &
            energy_diffusion,impact_ionization,clean_param,rhoe2T,T2rhoe, &
            rhov12v1,v12rhov1,clean_param_after_regrid
use advec, only: interface_vels_allspec,set_global_boundaries_allspec
use timeutils, only: dateinc
use io_nompi, only: interp_file2subgrid,plasma_output_nompi
use potential_nompi, only: set_fields_test,velocities_nompi
use geomagnetic, only: geog2geomag,ECEFspher2ENU
use interpolation, only: interp3,interp2
use calculus, only: grad3D2,grad3D3
use gemini_work_def, only: gemini_work

implicit none (type, external)
private
public :: c_params, gemini_alloc, gemini_dealloc, init_precipinput_in, msisinit_in, &
            set_start_values_auxtimevars, set_start_timefromcfg, set_start_values_auxvars, init_neutralBG_in, &
            set_update_cadence, neutral_atmos_winds, get_solar_indices, &
            v12rhov1_in, T2rhoe_in, interface_vels_allspec_in, sweep3_allparams_in, &
            sweep1_allparams_in, sweep2_allparams_in, &
            rhov12v1_in, VNRicht_artvisc_in, compression_in, rhoe2T_in, clean_param_in, &
            energy_diffusion_in, source_loss_allparams_in, &
            dateinc_in, get_subgrid_size,get_fullgrid_size,get_config_vars, get_species_size, fluidvar_pointers, &
            fluidauxvar_pointers, electrovar_pointers, gemini_work, &
            interp_file2subgrid_in,grid_from_extents_in,read_fullsize_gridcenter_in, &
            gemini_work_alloc, gemini_work_dealloc, gemini_cfg_alloc, cli_in, read_config_in, gemini_cfg_dealloc, &
            grid_size_in, gemini_double_alloc, gemini_double_dealloc, gemini_grid_alloc, gemini_grid_dealloc, &
            gemini_grid_generate, setv2v3, v2grid, v3grid, maxcfl_in, plasma_output_nompi_in, set_global_boundaries_allspec_in, &
            get_fullgrid_lims_in,checkE1,get_cfg_timevars,electrodynamics_test,forceZOH_all, permute_fluidvars, &
            ipermute_fluidvars, tag4refine_location, tag4refine_vperp, clean_param_after_regrid_in, get_locationsi_in, &
            set_datainow_in, get_datainow_ptr_in, swap_statevars, interp3_in, interp2_in, tag4refine_diff, &
            tag4refine_grad, tag4coarsening_diff


real(wp), protected :: v2grid,v3grid

!> type for passing C-like parameters between program units
type, bind(C) :: c_params
  !! this MUST match gemini3d.h and libgemini.f90 exactly including order
  logical(c_bool) :: fortran_nml, fortran_cli, debug, dryrun
  character(kind=c_char) :: out_dir(1000)
  !! .ini [base]
  integer(c_int) :: ymd(3)
  real(kind=c_float) :: UTsec0, tdur, dtout, activ(3), tcfl, Teinf
  !! .ini
end type c_params


interface checkarray
  procedure :: checkarray3D,checkarray4D
end interface checkarray


contains
  !> interface subroutine from which we can read in ONLY the grid sizes
  subroutine grid_size_in(cfg)
    type(gemini_cfg), intent(in) :: cfg

    call grid_size(cfg%indatsize)
  end subroutine grid_size_in


  !> interface subroutine to handle command line inputs or otherwise setup variables that would be specified
  !    from the command line.
  subroutine cli_in(p,lid2in,lid3in,cfg)
    type(c_params), intent(in) :: p
    integer, intent(inout) :: lid2in,lid3in
    type(gemini_cfg), intent(inout) :: cfg
    character(size(p%out_dir)) :: buf
    integer :: i

    if(p%fortran_cli) then
      call cli(cfg, lid2in, lid3in, debug)
    else
      buf = "" !< ensure buf has no garbage characters

      do i = 1, len(buf)
        if (p%out_dir(i) == c_null_char) exit
        buf(i:i) = p%out_dir(i)
      enddo
      cfg%outdir = expanduser(buf)

      cfg%dryrun = p%dryrun
      debug = p%debug
    endif
  end subroutine cli_in


  !> interface layer to find and read in the config file (we assume struct has already been allocated
  subroutine read_config_in(p,cfg)
    type(c_params), intent(in) :: p
    type(gemini_cfg), intent(inout) :: cfg

    !> read the config input file, if not passed .ini info from C++ frontend
    if(p%fortran_nml) then
      call find_config(cfg)
      call read_configfile(cfg, verbose=.false.)
    endif

    !> at this point we can check the input files and make sure we have a well-formed simulation setup
    call check_input_files(cfg)
  end subroutine read_config_in


  !> return some data from cfg that is needed in the main program
  subroutine get_config_vars(cfg,flagneuBG,flagdneu,dtneuBG,dtneu)
    type(gemini_cfg), intent(in) :: cfg
    logical, intent(inout) :: flagneuBG
    integer, intent(inout) :: flagdneu
    real(wp), intent(inout) :: dtneuBG,dtneu

    flagneuBG=cfg%flagneuBG; flagdneu=cfg%flagdneu;
    dtneuBG=cfg%dtneuBG; dtneu=cfg%dtneu;
  end subroutine get_config_vars


  !> returns the subgrid sizes *** stored in the grid module ***
  subroutine get_subgrid_size(lx1out,lx2out,lx3out)
    integer, intent(inout) :: lx1out,lx2out,lx3out

    lx1out=lx1; lx2out=lx2; lx3out=lx3;
  end subroutine get_subgrid_size


  !> return full grid extents *** stored in the grid module ***
  subroutine get_fullgrid_size(lx1out,lx2allout,lx3allout)
    integer, intent(inout) :: lx1out,lx2allout,lx3allout

    lx1out=lx1; lx2allout=lx2all; lx3allout=lx3all;
  end subroutine get_fullgrid_size


  !> return number of species *** from phys_consts module ***
  subroutine get_species_size(lspout)
    integer, intent(inout) :: lspout

    lspout=lsp
  end subroutine get_species_size


  !> return the limits of the grid to caller
  subroutine get_fullgrid_lims_in(x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax)
    real(wp), intent(inout) :: x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax

    call get_fullgrid_lims(x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax)
  end subroutine get_fullgrid_lims_in


  !> allocate space for config struct, and return a pointer
  function gemini_cfg_alloc() result(cfg)
     type(gemini_cfg), pointer :: cfg

     allocate(cfg)
  end function gemini_cfg_alloc


  !> deallocate config struct
  subroutine gemini_cfg_dealloc(cfg)
    type(gemini_cfg), pointer, intent(inout) :: cfg

    if (associated(cfg)) then
      deallocate(cfg)
      cfg=>null()
    end if
  end subroutine gemini_cfg_dealloc


  !> allocate struct for internal variables
  function gemini_work_alloc(cfg) result(intvars)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), pointer :: intvars

    !> none of this can be done unless the size variables in the grid module are set
    if (lx1<=0 .or. lx2<=0 .or. lx3<=0 .or. lx2all<=0 .or. lx3all<=0) then
      print*,  '  Malformed size from grid module:  ',lx1,lx2,lx3,lx2all,lx3all
      error stop
    end if

    allocate(intvars)

    !> neutral variables (never need to be haloed, etc.)
    allocate(intvars%atmos)
    call neutral_info_alloc(intvars%atmos)

     !> space for integrated volume emission rates (lx2,lx3,lwave)
    if (cfg%flagglow /= 0) then
      allocate(intvars%iver(lx2,lx3,lwave))
      intvars%iver = 0
    end if

    !> allocate space for some arrays needed for fluid solves, note that these arrays are not haloed; they
    !    are computed from haloed vs1,2,3 arrays
    allocate(intvars%vs1i(1:lx1+1,1:lx2,1:lx3,1:lsp))
    allocate(intvars%vs2i(1:lx1,1:lx2+1,1:lx3,1:lsp))
    allocate(intvars%vs3i(1:lx1,1:lx2,1:lx3+1,1:lsp))
    allocate(intvars%Q(1:lx1,1:lx2,1:lx3,1:lsp))
    intvars%vs1i=0._wp
    intvars%vs2i=0._wp
    intvars%vs3i=0._wp
    intvars%Q=0._wp

    allocate(intvars%eprecip)
    allocate(intvars%efield)
    ! fields of intvars%atmos are allocated in neutral:neutral_info_alloc()

    ! Here the user needs to allocate any custom variables they want to pass around and/or output
    allocate(intvars%sigP(1:lx1,1:lx2,1:lx3))
    allocate(intvars%sigH, mold=intvars%sigP)
  end function gemini_work_alloc


  !> deallocate struct for internal variables
  subroutine gemini_work_dealloc(cfg,intvars)
    type(gemini_cfg), intent(in) :: cfg
    type(gemini_work), pointer, intent(inout) :: intvars

    !> neutral variables (never need to be haloed, etc.)
    print*, 'Deallocating atmospheric state variables used in GEMINI...'
    call neutral_info_dealloc(intvars%atmos)
    deallocate(intvars%atmos)

    !> space for integrated volume emission rates (lx2,lx3,lwave)
    if (cfg%flagglow /= 0) then
      print*, 'Deallocating glow data:  '
      deallocate(intvars%iver)
    end if

    !> allocate space for some arrays needed for fluid solves, note that these arrays are not haloed; they
    !    are computed from haloed vs1,2,3 arrays
    print*, 'Deallocating internal variables for GEMINI...'
    deallocate(intvars%vs1i)
    deallocate(intvars%vs2i)
    deallocate(intvars%vs3i)
    deallocate(intvars%Q)

    if (associated(intvars%eprecip)) deallocate(intvars%eprecip)
    if (associated(intvars%efield)) deallocate(intvars%efield)
    !call clear_dneu(intvars%atmosperturb)    ! requies mpi so omitted here?

    if (associated(intvars%Phiall)) deallocate(intvars%Phiall)

    ! Here the user *must* deallocate their custom vars
    deallocate(intvars%sigP,intvars%sigH)

    deallocate(intvars)
  end subroutine gemini_work_dealloc


  !> allocate space for gemini state variables, bind pointers to blocks of memory; this is primarily meant
  !    to be called from a fortran main program and simply encapsulates a set of calls to other elementary
  !    allocation procedures which could alternatively be directly called from the main GEMINI app
  subroutine gemini_alloc(cfg,fluidvars,fluidauxvars,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars
    type(gemini_work), pointer, intent(inout) :: intvars

    !> allocate floating point arrays
    call gemini_double_alloc(fluidvars,fluidauxvars,electrovars)

    !> internal work struct
    intvars=>gemini_work_alloc(cfg)
  end subroutine gemini_alloc


  !> Fortran calls to allocate floating point arrays (should only be used from fortran)
  subroutine gemini_double_alloc(fluidvars,fluidauxvars,electrovars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars

    !> one contiguous block for overall simulation data
    allocate(fluidvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,5*lsp))
    !> fluid momentum and energy density variables
    allocate(fluidauxvars(-1:lx1+2,-1:lx2+2,-1:lx3+2,2*lsp+9))
    !> electrodynamic state variables (lx1,lx2,lx3)
    allocate(electrovars(-1:lx1+2,-1:lx2+2,-1:lx3+2,7))

    !> this is a safety bit of code to make sure everything starts to zero; apparently some things are not getting initialized in some cases
    fluidvars=0._wp
    fluidauxvars=0._wp
    electrovars=0._wp
  end subroutine gemini_double_alloc


  subroutine gemini_double_dealloc(fluidvars,fluidauxvars,electrovars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: electrovars

    !> ifort generates a runtime error for this if called from C; I guess memory management needs to be done on the C-side of things
    deallocate(fluidvars)
    deallocate(fluidauxvars)
    deallocate(electrovars)
  end subroutine gemini_double_dealloc


  !> sequence of calls to allocate space for grid object and internal variables (analogous to gemini_work_alloc)
  subroutine gemini_grid_alloc(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x,xtype,xC)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout), pointer :: x
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    integer :: ix2,ix3
    real(wp), dimension(:), allocatable :: x2,x3
    real(wp), dimension(:), pointer :: x1
    real(wp) :: glonctr,glatctr

    ! retrieve data from grid module
    call get_gridcenter(glonctr,glatctr)

    ! error checking
    if (glatctr<-90._wp .or. glatctr>90._wp) then
      error stop ' grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter to assign &
                   module variables glonctr,glatctr'
    end if

    ! create temp space
    allocate(x2(lx2wg),x3(lx3wg))

    ! make uniformly spaced coordinate arrays; unless a x1 array was already provided by user
    call get_x1coords(x1)
    !x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)]

    ! allocate mesh class and create a C pointer to it
    call meshobj_alloc(x1,x2,x3,x,xtype,xC)

    ! allocate grid internal data arrays/structs
    call grid_internaldata_alloc(x1,x2,x3,x2,x3,glonctr,glatctr,x)

    ! get rid of temp. arrays
    deallocate(x2,x3)
  end subroutine gemini_grid_alloc


  !> subroutine to force generate of grid internal objects (grid must already be allocated)
  subroutine gemini_grid_generate(x)
    class(curvmesh), intent(inout) :: x

    call grid_internaldata_generate(x)
  end subroutine gemini_grid_generate


  !> deallocate grid data
  subroutine gemini_grid_dealloc(x,xtype,xC)
    class(curvmesh), pointer, intent(inout) :: x
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC

    call grid_internaldata_ungenerate(x)    ! this both ungenerates and also deallocates data stored in grid object
    call meshobj_dealloc(x,xtype,xC)
  end subroutine gemini_grid_dealloc


  !> force a value for the lagrangian grid drift
  subroutine setv2v3(v2gridin,v3gridin)
    real(wp), intent(in) :: v2gridin,v3gridin

    v2grid=v2gridin
    v3grid=v3gridin
  end subroutine setv2v3


  !> take a block of memory and assign pointers to various pieces representing different fluid, etc. state variables
  !!   This will be called any time a gemini library procedures needs to access individual state variables.
  subroutine fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: ns
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: vs1,vs2,vs3
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: Ts

    if (.not. associated(fluidvars)) error stop ' Attempting to bind fluid state vars to unassociated memory!'

    !> main state variables for gemini (lx1+4,lx2+4,lx3+4,lsp)
    ns=>fluidvars(:,:,:,1:lsp)
    vs1=>fluidvars(:,:,:,lsp+1:2*lsp)
    vs2=>fluidvars(:,:,:,2*lsp+1:3*lsp)
    vs3=>fluidvars(:,:,:,3*lsp+1:4*lsp)
    Ts=>fluidvars(:,:,:,4*lsp+1:5*lsp)
  end subroutine


  !> bind pointers for auxiliary fluid variables to a contiguous block of memory
  subroutine fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer, intent(inout) :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    if (.not. associated(fluidauxvars)) error stop ' Attempting to bind aux fluid state vars to unassociated memory!'

    !> pointers to aliased state variables
    rhovs1=>fluidauxvars(:,:,:,1:lsp)
    rhoes=>fluidauxvars(:,:,:,lsp+1:2*lsp)

    !> MHD-like state variables used in some calculations (lx1+4,lx2+4,lx3+4,lsp)
    rhov2=>fluidauxvars(:,:,:,2*lsp+1)
    rhov3=>fluidauxvars(:,:,:,2*lsp+2)
    B1=>fluidauxvars(:,:,:,2*lsp+3)
    B2=>fluidauxvars(:,:,:,2*lsp+4)
    B3=>fluidauxvars(:,:,:,2*lsp+5)
    v1=>fluidauxvars(:,:,:,2*lsp+6)
    v2=>fluidauxvars(:,:,:,2*lsp+7)
    v3=>fluidauxvars(:,:,:,2*lsp+8)
    rhom=>fluidauxvars(:,:,:,2*lsp+9)
  end subroutine fluidauxvar_pointers


  !> bind pointers for electomagnetic state variables to a contiguous block of memory
  subroutine electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    real(wp), dimension(:,:,:), pointer, intent(inout) :: E1,E2,E3,J1,J2,J3,Phi

    if (.not. associated(electrovars)) error stop ' Attempting to bind electro state vars to unassociated memory!'

    !> electric fields, potential, and current density
    E1=>electrovars(:,:,:,1)
    E2=>electrovars(:,:,:,2)
    E3=>electrovars(:,:,:,3)
    J1=>electrovars(:,:,:,4)
    J2=>electrovars(:,:,:,5)
    J3=>electrovars(:,:,:,6)
    Phi=>electrovars(:,:,:,7)
  end subroutine electrovar_pointers


  !> deallocate state variables include double precision data arrays; only works for all compilers from fortran main programs
  !    This is a wrapper that successively calls all deallocations and it mean to be used from a fortran app.
  subroutine gemini_dealloc(cfg,fluidvars,fluidauxvars,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars, fluidauxvars, electrovars
    type(gemini_work), pointer, intent(inout) :: intvars

    !> ifort generates a runtime error for this if called from C; I guess memory management needs to be done on the C-side of things
    call gemini_double_dealloc(fluidvars,fluidauxvars,electrovars)

    !call gemini_dealloc_nodouble(cfg,intvars)
    call gemini_work_dealloc(cfg,intvars)
  end subroutine


  !> Have each worker read the entire input file and then interpolate it onto its own subgrid
  subroutine interp_file2subgrid_in(cfg,x,fluidvars,electrovars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars,electrovars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    logical :: errflag=.false.
    real(wp) :: nlower,nupper,vlower,vupper,Tlower,Tupper

    print*, 'Initiating tiled/interpolate input...'
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call interp_file2subgrid(cfg%indatsize,cfg%indatfile,cfg%indatgrid,x%x1,x%x2,x%x3,ns,vs1,Ts,Phi)

    ! it's important to note that the file input *does not specify* vs2,vs3 so we need to set some initial
    !   values here, otherwise we may just get existing garbage in memory.
    vs2=0._wp; vs3=0._wp;

    ! the input code also does not assign electrodynamic variables
    E1=0._wp; E2=0._wp; E3=0._wp; J1=0._wp; J2=0._wp; J3=0._wp;

    ! Clean the input data since it has been interpolated and there is a chance we get unacceptable values
    !   This appears to sometimes be needed with dipole grid simulations.
    call clean_param(x,1,ns)
    call clean_param(x,2,vs1)
    call clean_param(x,3,Ts)

    ! this is a good time to do some error checking since each patch only calls this code once per simulation
    nlower=0; nupper=1e14;
    vlower=-1e4; vupper=1e4;
    Tlower=0; Tupper=20000;
    print*, 'In fortran file i/o...'
    errflag=errflag .or. checkarray(ns,nlower,nupper,'>>> Full density corrupted:  ',0)
    errflag=errflag .or. checkarray(vs1,vlower,vupper,'>>> Full velocity corrupted:  ',0)
    errflag=errflag .or. checkarray(Ts,Tlower,Tupper,'>>> Full temperature corrupted:  ',0)
    if (errflag) error stop
    print*, 'Exiting fortran file i/o...'
  end subroutine interp_file2subgrid_in


  !> Basic utility to have each worker dump state variable contents to a file
  subroutine plasma_output_nompi_in(cfg,ymd,UTsec,fluidvars,electrovars,identifier,x1lims,x2lims,x3lims)
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars,electrovars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    integer, intent(in) :: identifier
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call plasma_output_nompi(cfg%outdir,cfg%flagoutput,ymd,UTsec,ns, &
                               vs1,vs2,vs3,Ts,Phi,J1,J2,J3, &
                               identifier,x1lims,x2lims,x3lims)
  end subroutine plasma_output_nompi_in


  !> interface for pulling grid center coordinates from the input file
  subroutine read_fullsize_gridcenter_in(cfg)
    type(gemini_cfg), intent(in) :: cfg

    call read_size_gridcenter(cfg%indatsize,cfg%indatgrid)
  end subroutine read_fullsize_gridcenter_in


  !> A somewhat superfluous wrapper for grid generation from known extents, included here to keep with the
  !    pattern of only having applications access procedures from this module and no others.  In the case
  !    of a Cartesian grid we need to somehow set the center point so that there are glat/glon associated
  !    with each location of the mesh.  This will be obtained from the input cfg file in cases where we are
  !    generating a grid from extents.  Note that this is distinctly different from the situation where we
  !    are using read_grid() to input a grid from a file - in that case the parameters glonctr and glatctr
  !    are expected to be kept within that file.
  subroutine grid_from_extents_in(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout) :: x

    call grid_from_extents(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
  end subroutine grid_from_extents_in


  !> assign initial values on some auxiliary time variables
  subroutine set_start_values_auxtimevars(it,t,tout,tglowout,tneuBG)
    integer, intent(inout) :: it
    real(wp), intent(inout) :: t,tout,tglowout,tneuBG

    !> Initialize some variables need for time stepping and output
!    it = 1; t = 0; tout = t; tglowout = t; tneuBG=t
    it = 1; tout = t; tglowout = t; tneuBG=t
  end subroutine set_start_values_auxtimevars


  ! pull relevant time variables from the cfg structure
  subroutine get_cfg_timevars(cfg,tmilestone,flagneuBG,dtneuBG,flagdneu,flagoutput)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(inout) :: tmilestone
    logical, intent(inout) :: flagneuBG
    real(wp), intent(inout) :: dtneuBG
    integer, intent(inout) :: flagdneu
    integer, intent(inout) :: flagoutput

    tmilestone=0._wp   ! make sure first output is a milestone
    flagneuBG=cfg%flagneuBG
    dtneuBG=cfg%dtneuBG
    flagdneu=cfg%flagdneu
    flagoutput=cfg%flagoutput
  end subroutine get_cfg_timevars


  !> Force time values to a specific date/time and set duration based on cfg argument fields
  subroutine set_start_timefromcfg(cfg,ymd,UTsec,tdur)
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec
    real(wp), intent(inout) :: tdur

    UTsec = cfg%UTsec0
    ymd = cfg%ymd0
    tdur = cfg%tdur
  end subroutine set_start_timefromcfg


  !> set start values for some variables not specified by the input files.
  !    some care is required here because the state variable pointers are mapped;
  !    however, note that the lbound and ubound have not been set since arrays are not passed through as dummy args
  !    with specific ubound so that we need to use intrinsic calls to make sure we fill computational cells (not ghost)
  subroutine set_start_values_auxvars(x,fluidauxvars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
    !THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
    !FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INERTIAL CAPACITANCE THERE
    !WILL BE A FINITE AMOUNT OF TIME FOR THE FLOWS TO 'START UP', BUT THIS SHOULDN'T
    !BE TOO MUCH OF AN ISSUE.  WE ALSO NEED TO SET THE BACKGROUND MAGNETIC FIELD STATE
    !VARIABLE HERE TO WHATEVER IS SPECIFIED IN THE GRID STRUCTURE (THESE MUST BE CONSISTENT)
    rhov2 = 0; rhov3 = 0; v2 = 0; v3 = 0; B2 = 0; B3 = 0;
    call set_magfield(x,fluidauxvars)
  end subroutine set_start_values_auxvars


  subroutine set_magfield(x,fluidauxvars)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    integer :: ix1min,ix1max,ix2min,ix2max,ix3min,ix3max
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    ix1min=lbound(B1,1)+2
    ix1max=ubound(B1,1)-2
    ix2min=lbound(B1,2)+2
    ix2max=ubound(B1,2)-2
    ix3min=lbound(B1,3)+2
    ix3max=ubound(B1,3)-2
    B1(ix1min:ix1max,ix2min:ix2max,ix3min:ix3max) = x%Bmag(1:lx1,1:lx2,1:lx3)
    !! this assumes that the grid is defined s.t. the x1 direction corresponds
    !! to the magnetic field direction (hence zero B2 and B3).
  end subroutine set_magfield


  !> Wrapper for initialization of electron precipitation data
  subroutine init_precipinput_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call init_precipinput(dt, cfg,ymd,UTsec,x,intvars%eprecip)
  end subroutine init_precipinput_in


  !> initialization procedure needed for MSIS 2.0
  subroutine msisinit_in(cfg)
    type(gemini_cfg), intent(in) :: cfg
    logical :: exists

    character(len=11) :: msis2_param_file

    select case (cfg%msis_version)
    case(0)
      !! MSISE00
      return
    case(20)
      msis2_param_file = "msis20.parm"
    case(21)
      msis2_param_file = "msis21.parm"
    case default
      !! new or unknown version of MSIS, default MSIS 2.x parameter file
      msis2_param_file = ""
    end select

    if(len_trim(msis2_param_file) > 0) then
      inquire(file=msis2_param_file, exist=exists)
      if(.not.exists) error stop 'could not find MSIS 2 parameter file ' // msis2_param_file // &
        ' this file must be in the same directory as gemini.bin, and run from that directory. ' // &
        'This limitation comes from how MSIS 2.x is coded internally.'
      call msisinit(parmfile=msis2_param_file)
    else
      call msisinit()
    end if
  end subroutine msisinit_in


  !> call to initialize the neutral background information
  subroutine init_neutralBG_in(cfg,x,dt,t,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(inout) :: x    ! so neutral module can deallocate unit vectors once used...
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call init_neutralBG(cfg,ymd,UTsec,x,v2grid,v3grid,intvars%atmos)
  end subroutine init_neutralBG_in


  !> set update cadence for printing out diagnostic information during simulation
  subroutine set_update_cadence(iupdate)
    integer, intent(inout) :: iupdate

    !> control update rate from excessive console printing
    !! considering small vs. large simulations
    !! these are arbitrary levels, so feel free to finesse
    if (lx1*lx2*lx3 < 20000) then
      iupdate = 50
    elseif (lx1*lx2*lx3 < 100000) then
      iupdate = 10
    else
      iupdate = 1
    endif
  end subroutine set_update_cadence


  !> compute background neutral density, temperature, and wind
  subroutine neutral_atmos_winds(cfg,x,ymd,UTsec,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(gemini_work), intent(inout) :: intvars

    call neutral_atmos(ymd,UTsec,x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3), &
                         cfg%activ,cfg%msis_version, atmos=intvars%atmos)
    call neutral_winds(ymd, UTsec, Ap=cfg%activ(3), x=x, atmos=intvars%atmos)
  end subroutine neutral_atmos_winds


  !> get solar indices from cfg struct
  subroutine get_solar_indices(cfg,f107,f107a)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(inout) :: f107,f107a

    f107=cfg%activ(2)
    f107a=cfg%activ(1)
  end subroutine get_solar_indices


  !> convert velocity to momentum density
  subroutine v12rhov1_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call v12rhov1(ns,vs1,rhovs1)
  end subroutine v12rhov1_in


  !> convert temperature to specific internal energy density
  subroutine T2rhoe_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call T2rhoe(ns,Ts,rhoes)
  end subroutine T2rhoe_in


  !> compute interface velocities once haloing has been done
  subroutine interface_vels_allspec_in(fluidvars,intvars,lsp)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    type(gemini_work), intent(inout) :: intvars
    integer, intent(in) :: lsp
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call interface_vels_allspec(vs1,vs2,vs3,intvars%vs1i,intvars%vs2i,intvars%vs3i,lsp)    ! needs to happen regardless of ions v. electron due to energy eqn.
  end subroutine interface_vels_allspec_in


  !> enforce global boundary conditions
  subroutine set_global_boundaries_allspec_in(x,fluidvars,fluidauxvars,intvars,lsp)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    integer, intent(in) :: lsp

    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    ! bind pointers
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    ! fix global boundaries, as needed
    call set_global_boundaries_allspec(x%flagper,ns,rhovs1,vs1,vs2,vs3,rhoes,intvars%vs1i,lsp,x)
  end subroutine set_global_boundaries_allspec_in


  !> functions for sweeping advection
  subroutine sweep3_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep3_allparams(dt,x,intvars%vs3i,ns,rhovs1,rhoes)
  end subroutine sweep3_allparams_in


  subroutine sweep1_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep1_allparams(dt,x,intvars%vs1i,ns,rhovs1,rhoes)
  end subroutine sweep1_allparams_in


  subroutine sweep2_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(inout) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call sweep2_allparams(dt,x,intvars%vs2i,ns,rhovs1,rhoes)
  end subroutine sweep2_allparams_in


  !> conversion of momentum density to velocity
  subroutine rhov12v1_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call rhov12v1(ns,rhovs1,vs1)
  end subroutine rhov12v1_in


  !> compute artifical viscosity
  subroutine VNRicht_artvisc_in(fluidvars,intvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    type(gemini_work), intent(inout) :: intvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call VNRicht_artvisc(ns,vs1,intvars%Q)
  end subroutine VNRicht_artvisc_in


  !> compression substep for fluid solve
  subroutine compression_in(fluidvars,fluidauxvars,intvars,x,dt)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    type(gemini_work), intent(in) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    call compression(dt,x,vs1,vs2,vs3,intvars%Q,rhoes)   ! this applies compression substep
  end subroutine compression_in


  !> convert specific internal energy density into temperature
  subroutine rhoe2T_in(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call rhoe2T(ns,rhoes,Ts)
  end subroutine rhoe2T_in


  !> deal with null cell solutions
  subroutine clean_param_in(iparm,x,fluidvars)
    integer, intent(in) :: iparm
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: parm
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    select case (iparm)
      case (1)
        parm=>ns
      case (2)
        parm=>vs1
      case (3)
        parm=>Ts
      case default
        error stop '  libgemini:clean_params_C(); invalid parameter selected'
    end select
    call clean_param(x,iparm,parm)
  end subroutine clean_param_in


  !> deal with null cell solutions
  subroutine clean_param_after_regrid_in(iparm,x,fluidvars,intvars)
    integer, intent(in) :: iparm
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    type(gemini_work), intent(in) :: intvars
    real(wp), dimension(:,:,:,:), pointer :: parm
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    select case (iparm)
      case (1)
        parm=>ns
      case (2)
        parm=>vs1
      case (3)
        parm=>Ts
      case default
        error stop '  libgemini:clean_params_C(); invalid parameter selected'
    end select
    call clean_param_after_regrid(x,iparm,parm,intvars%atmos%Tn)
  end subroutine clean_param_after_regrid_in


  !> diffusion of energy
  subroutine energy_diffusion_in(cfg,x,fluidvars,electrovars,intvars,dt)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(in) :: intvars
    real(wp), intent(in) :: dt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call energy_diffusion(dt,x,ns,Ts,J1,intvars%atmos%nn,intvars%atmos%Tn,cfg%diffsolvetype,cfg%Teinf)
  end subroutine energy_diffusion_in


  !> source/loss numerical solutions
  subroutine source_loss_allparams_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt,t,ymd, &
                                        UTsec,f107a,f107,first,gavg,Tninf)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    type(gemini_work), intent(in) :: intvars
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: f107a,f107
    logical, intent(in) :: first
    real(wp), intent(in) :: gavg,Tninf
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    call source_loss_allparams(dt,t,cfg,ymd,UTsec,x,E1,E2,E3,intvars%Q,f107a,f107,intvars%atmos%nn, &
                                     intvars%atmos%vn1,intvars%atmos%vn2,intvars%atmos%vn3, &
                                     intvars%atmos%Tn,first,ns,rhovs1,rhoes,vs1,vs2,vs3,Ts, &
                                     intvars%iver,gavg,Tninf,intvars%eprecip)
  end subroutine source_loss_allparams_in


  !> print out min/max values for variables
  subroutine checkE1(fluidvars,fluidauxvars,electrovars,locID)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: electrovars
    integer, intent(in) :: locID
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:),pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: nlower,nupper,vlower,vupper,Tlower,Tupper
    real(wp) :: vplower,vpupper
    real(wp) :: Eparlower,Eparupper,Elower,Eupper,Jlower,Jupper,Philower,Phiupper
    real(wp) :: rhovlower,rhovupper,rhoeslower,rhoesupper,Blower,Bupper
    logical :: errflag=.false.
    integer :: funit

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check main variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! force a check on the interior cells of the domain
    nlower=0; nupper=1e14;
    vlower=-1e4; vupper=1e4;
    vplower=-1e4; vpupper=1e4;
    Tlower=0; Tupper=100000;

    errflag=errflag .or. checkarray(ns(3:lx1+2,3:lx2+2,3:lx3+2,:),nlower,nupper, &
                                     '>>> Interior density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(3:lx1+2,3:lx2+2,3:lx3+2,:),vlower,vupper, &
                                     '>>> Interior velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(3:lx1+2,3:lx2+2,3:lx3+2,:),Tlower,Tupper, &
                                     '>>> Interior temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:),vplower,vpupper, &
                                     '>>> Interior v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(3:lx1+2,3:lx2+2,3:lx3+2,:),vplower,vpupper, &
                                     '>>> Interior v3 data corrupted:  ',locID)

    ! now check the bottom ghost cells
    errflag=errflag .or. checkarray(ns(1:2,:,:,:),nlower,nupper, &
                                     '>>> Bottom density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(1:2,:,:,:),vlower,vupper, &
                                     '>>> Bottom velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(1:2,:,:,:),Tlower,Tupper, &
                                     '>>> Bottom temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(1:2,:,:,:),vplower,vpupper, &
                                     '>>> Bottom v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(1:2,:,:,:),vplower,vpupper, &
                                     '>>> Bottom v3 data corrupted:  ',locID)

    ! check top
    errflag=errflag .or. checkarray(ns(lx1+3:lx1+4,:,:,:),nlower,nupper, &
                                     '>>> Top density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(lx1+3:lx1+4,:,:,:),vlower,vupper, &
                                     '>>> Top velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(lx1+3:lx1+4,:,:,:),Tlower,Tupper, &
                                     '>>> Top temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(lx1+3:lx1+4,:,:,:),vplower,vpupper, &
                                     '>>> Top v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(lx1+3:lx1+4,:,:,:),vplower,vpupper, &
                                     '>>> Top v3 data corrupted:  ',locID)

    ! check left
    errflag=errflag .or. checkarray(ns(:,1:2,:,:),nlower,nupper, &
                                     '>>> Left density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,1:2,:,:),vlower,vupper, &
                                     '>>> Left velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,1:2,:,:),Tlower,Tupper, &
                                     '>>> Left temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,1:2,:,:),vplower,vpupper, &
                                     '>>> Left v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,1:2,:,:),vplower,vpupper, &
                                     '>>> Left v3 data corrupted:  ',locID)

    ! check right
    errflag=errflag .or. checkarray(ns(:,lx2+3:lx2+4,:,:),nlower,nupper, &
                                     '>>> Right density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,lx2+3:lx2+4,:,:),vlower,vupper, &
                                     '>>> Right velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,lx2+3:lx2+4,:,:),Tlower,Tupper, &
                                     '>>> Right temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,lx2+3:lx2+4,:,:),vplower,vpupper, &
                                     '>>> Right v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,lx2+3:lx2+4,:,:),vplower,vpupper, &
                                     '>>> Right v3 data corrupted:  ',locID)

    ! check bwd
    errflag=errflag .or. checkarray(ns(:,:,1:2,:),nlower,nupper, &
                                     '>>> Bwd density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,:,1:2,:),vlower,vupper, &
                                     '>>> Bwd velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,:,1:2,:),Tlower,Tupper, &
                                     '>>> Bwd temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,:,1:2,:),vplower,vpupper, &
                                     '>>> Bwd v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,:,1:2,:),vplower,vpupper, &
                                     '>>> Bwd v3 data corrupted:  ',locID)

    ! check fwd
    errflag=errflag .or. checkarray(ns(:,:,lx3+3:lx3+4,:),nlower,nupper, &
                                     '>>> Fwd density data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs1(:,:,lx3+3:lx3+4,:),vlower,vupper, &
                                     '>>> Fwd velocity data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Ts(:,:,lx3+3:lx3+4,:),Tlower,Tupper, &
                                     '>>> Fwd temperature data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs2(:,:,lx3+3:lx3+4,:),vplower,vpupper, &
                                     '>>> Fwd v2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(vs3(:,:,lx3+3:lx3+4,:),vplower,vpupper, &
                                     '>>> Fwd v3 data corrupted:  ',locID)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check electro variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Eparlower=-1e-4; Eparupper=1e-4;
    Elower=-0.5; Eupper=0.5;
    Jlower=-1e-3; Jupper=1e-3;
    Philower=-1e5; Phiupper=1e5;

    errflag=errflag .or. checkarray(E1(3:lx1+2,3:lx2+2,3:lx3+2),Eparlower,Eparupper, &
                                     '>>> Interior E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(3:lx1+2,3:lx2+2,3:lx3+2),Elower,Eupper, &
                                     '>>> Interior E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(3:lx1+2,3:lx2+2,3:lx3+2),Elower,Eupper, &
                                     '>>> Interior E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(3:lx1+2,3:lx2+2,3:lx3+2),Jlower,Jupper, &
                                     '>>> Interior J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(3:lx1+2,3:lx2+2,3:lx3+2),Jlower,Jupper, &
                                     '>>> Interior J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(3:lx1+2,3:lx2+2,3:lx3+2),Jlower,Jupper, &
                                     '>>> Interior J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(3:lx1+2,3:lx2+2,3:lx3+2),Philower,Phiupper, &
                                     '>>> Interior Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(1:2,:,:),Eparlower,Eparupper, &
                                     '>>> Bottom E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(1:2,:,:),Elower,Eupper, &
                                     '>>> Bottom E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(1:2,:,:),Elower,Eupper, &
                                     '>>> Bottom E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(1:2,:,:),Jlower,Jupper, &
                                     '>>> Bottom J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(1:2,:,:),Jlower,Jupper, &
                                     '>>> Bottom J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(1:2,:,:),Jlower,Jupper, &
                                     '>>> Bottom J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(1:2,:,:),Philower,Phiupper, &
                                     '>>> Bottom Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(lx1+3:lx1+4,:,:),Eparlower,Eparupper, &
                                     '>>> Top E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(lx1+3:lx1+4,:,:),Elower,Eupper, &
                                     '>>> Top E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(lx1+3:lx1+4,:,:),Elower,Eupper, &
                                     '>>> Top E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(lx1+3:lx1+4,:,:),Jlower,Jupper, &
                                     '>>> Top J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(lx1+3:lx1+4,:,:),Jlower,Jupper, &
                                     '>>> Top J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(lx1+3:lx1+4,:,:),Jlower,Jupper, &
                                     '>>> Top J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(lx1+3:lx1+4,:,:),Philower,Phiupper, &
                                     '>>> Top Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,1:2,:),Eparlower,Eparupper, &
                                     '>>> Left E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,1:2,:),Elower,Eupper, &
                                     '>>> Left E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,1:2,:),Elower,Eupper, &
                                     '>>> Left E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,1:2,:),Jlower,Jupper, &
                                     '>>> Left J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,1:2,:),Jlower,Jupper, &
                                     '>>> Left J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,1:2,:),Jlower,Jupper, &
                                     '>>> Left J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,1:2,:),Philower,Phiupper, &
                                     '>>> Left Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,lx2+3:lx2+4,:),Eparlower,Eparupper, &
                                     '>>> Right E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,lx2+3:lx2+4,:),Elower,Eupper, &
                                     '>>> Right E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,lx2+3:lx2+4,:),Elower,Eupper, &
                                     '>>> Right E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,lx2+3:lx2+4,:),Jlower,Jupper, &
                                     '>>> Right J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,lx2+3:lx2+4,:),Jlower,Jupper, &
                                     '>>> Right J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,lx2+3:lx2+4,:),Jlower,Jupper, &
                                     '>>> Right J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,lx2+3:lx2+4,:),Philower,Phiupper, &
                                     '>>> Right Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,:,1:2),Eparlower,Eparupper, &
                                     '>>> Bwd E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,:,1:2),Elower,Eupper, &
                                     '>>> Bwd E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,:,1:2),Elower,Eupper, &
                                     '>>> Bwd E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,:,1:2),Jlower,Jupper, &
                                     '>>> Bwd J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,:,1:2),Jlower,Jupper, &
                                     '>>> Bwd J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,:,1:2),Jlower,Jupper, &
                                     '>>> Bwd J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,:,1:2),Philower,Phiupper, &
                                     '>>> Bwd Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,:,lx3+3:lx3+4),Eparlower,Eparupper, &
                                     '>>> Fwd E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd Phi data corrupted:  ',locID)

    errflag=errflag .or. checkarray(E1(:,:,lx3+3:lx3+4),Eparlower,Eparupper, &
                                     '>>> Fwd E1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E2(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(E3(:,:,lx3+3:lx3+4),Elower,Eupper, &
                                     '>>> Fwd E3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J1(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J2(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J2 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(J3(:,:,lx3+3:lx3+4),Jlower,Jupper, &
                                     '>>> Fwd J3 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(Phi(:,:,lx3+3:lx3+4),Philower,Phiupper, &
                                     '>>> Fwd Phi data corrupted:  ',locID)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check aux variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rhoeslower=0; rhoesupper=1;
    rhovlower=-1; rhovupper=1;
    Blower=-100000e-9; Bupper=100000e-9;

    errflag=errflag .or. checkarray(rhovs1(3:lx1+2,3:lx2+2,3:lx3+2,:),rhovlower,rhovupper, &
                                     '>>> Interior rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(3:lx1+2,3:lx2+2,3:lx3+2,:),rhoeslower,rhoesupper, &
                                     '>>> Interior rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(3:lx1+2,3:lx2+2,3:lx3+2),Blower,Bupper, &
                                     '>>> Interior B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(1:2,:,:,:),rhovlower,rhovupper, &
                                     '>>> Bottom rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(1:2,:,:,:),rhoeslower,rhoesupper, &
                                     '>>> Bottom rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(1:2,:,:),Blower,Bupper, &
                                     '>>> Bottom B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(lx1+3:lx1+4,:,:,:),rhovlower,rhovupper, &
                                     '>>> Top rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(lx1+3:lx1+4,:,:,:),rhoeslower,rhoesupper, &
                                     '>>> Top rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(lx1+3:lx1+4,:,:),Blower,Bupper, &
                                     '>>> Top B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,1:2,:,:),rhovlower,rhovupper, &
                                     '>>> Left rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,1:2,:,:),rhoeslower,rhoesupper, &
                                     '>>> Left rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,1:2,:),Blower,Bupper, &
                                     '>>> Left B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,lx2+3:lx2+4,:,:),rhovlower,rhovupper, &
                                     '>>> Right rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,lx2+3:lx2+4,:,:),rhoeslower,rhoesupper, &
                                     '>>> Right rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,lx2+3:lx2+4,:),Blower,Bupper, &
                                     '>>> Right B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,:,1:2,:),rhovlower,rhovupper, &
                                     '>>> Bwd rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,:,1:2,:),rhoeslower,rhoesupper, &
                                     '>>> Bwd rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,:,1:2),Blower,Bupper, &
                                     '>>> Bwd B1 data corrupted:  ',locID)

    errflag=errflag .or. checkarray(rhovs1(:,:,lx3+3:lx3+4,:),rhovlower,rhovupper, &
                                     '>>> Fwd rhovs1 data corrupted:  ',locID)
    errflag=errflag .or. checkarray(rhoes(:,:,lx3+3:lx3+4,:),rhoeslower,rhoesupper, &
                                     '>>> Fwd rhoes data corrupted:  ',locID)
    errflag=errflag .or. checkarray(B1(:,:,lx3+3:lx3+4),Blower,Bupper, &
                                     '>>> Fwd B1 data corrupted:  ',locID)

    if (errflag) then
      open(newunit=funit,file='error.dat',status='replace',access='stream')
      write(funit) ns
      write(funit) vs1
      write(funit) Ts
      close(funit)

      error stop
    end if
  end subroutine checkE1


  !> For purposes of testing we just want to set some values for the electric fields and comput drifts.
  subroutine electrodynamics_test(cfg,x,fluidvars,fluidauxvars,electrovars,intvars)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work), intent(in) :: intvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    integer :: lx1,lx2,lx3,lsp
    real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav
    real(wp), dimension(:,:,:,:), allocatable :: muP,muH,nusn

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
    lx1=x%lx1; lx2=x%lx2; lx3=x%lx3; lsp=size(ns,4);

!    print*, '-----------------------------------------------------------------------'
!!    print*, '                      ',minval(Ts(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(Ts(3:lx1+2,3:lx2+2,3:lx3+2,:))
!!    print*, '                      ',minval(ns(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(ns(3:lx1+2,3:lx2+2,3:lx3+2,:))
!    print*, '                      ',minval(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:)), &
!                                       minval(vs3(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(vs3(3:lx1+2,3:lx2+2,3:lx3+2,:)), &
!                                       minval(vs1(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(vs1(3:lx1+2,3:lx2+2,3:lx3+2,:))
!!    print*, '                      ',maxloc(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:))
!!    print*, '                      ',shape(vs2),lbound(vs2)
!    print*, '                      ',minval(intvars%atmos%vn2),maxval(intvars%atmos%vn2), &
!                                       minval(intvars%atmos%vn3),maxval(intvars%atmos%vn3)
!    print*, '-----------------------------------------------------------------------'

    !call set_fields_test(x,E1,E2,E3)
    allocate(sig0(lx1,lx2,lx3),sigP(lx1,lx2,lx3),sigH(lx1,lx2,lx3),sigPgrav(lx1,lx2,lx3),sigHgrav(lx1,lx2,lx3))
    allocate(muP(lx1,lx2,lx3,lsp),muH(lx1,lx2,lx3,lsp),nusn(lx1,lx2,lx3,lsp))
    call conductivities(intvars%atmos%nn,intvars%atmos%Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
    call velocities_nompi(muP,muH,nusn,E2,E3,intvars%atmos%vn2,intvars%atmos%vn3,ns,Ts,x, &
                      cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
    deallocate(sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)

!    print*, '======================================================================='
!!    print*, '                      ',minval(Ts(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(Ts(3:lx1+2,3:lx2+2,3:lx3+2,:))
!!    print*, '                      ',minval(ns(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(ns(3:lx1+2,3:lx2+2,3:lx3+2,:))
!    print*, '                      ',minval(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:)), &
!                                       minval(vs3(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(vs3(3:lx1+2,3:lx2+2,3:lx3+2,:)), &
!                                       minval(vs1(3:lx1+2,3:lx2+2,3:lx3+2,:)),maxval(vs1(3:lx1+2,3:lx2+2,3:lx3+2,:))
!!    print*, '                      ',maxloc(vs2(3:lx1+2,3:lx2+2,3:lx3+2,:))
!    print*, '                      ',minval(intvars%atmos%vn2),maxval(intvars%atmos%vn2), &
!                                       minval(intvars%atmos%vn3),maxval(intvars%atmos%vn3)
!    print*, '======================================================================='
!!    error stop
  end subroutine electrodynamics_test


  !> utility procedure to check that array values are in a certain bounds
  logical function checkarray4D(array,lower,upper,errmsg,locID) result(errflag)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), intent(in) :: lower,upper
    character(*), intent(in) :: errmsg
    integer, intent(in) :: locID

    if (minval(array) < lower .or. maxval(array) > upper) then
      print*, locID,' ',errmsg,minval(array),maxval(array),minloc(array),maxloc(array)
      errflag=.true.
    else
      errflag=.false.
    end if
  end function checkarray4D


  logical function checkarray3D(array,lower,upper,errmsg,locID) result(errflag)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), intent(in) :: lower,upper
    character(*), intent(in) :: errmsg
    integer, intent(in) :: locID

    if (minval(array) < lower .or. maxval(array) > upper) then
      print*, locID,' ',errmsg,minval(array),maxval(array),minloc(array),maxloc(array)
      errflag=.true.
    else
      errflag=.false.
    end if
  end function checkarray3D


  !> return the maximum cfl over the grid
  subroutine maxcfl_in(fluidvars,x,dt,maxcfl)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars
    class(curvmesh), pointer :: x
    real(wp), intent(in) :: dt
    real(wp), intent(inout) :: maxcfl
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call cflcalc(Ts,vs1,vs2,vs3,x%dl1i,x%dl2i,x%dl3i,dt,maxcfl)
  end subroutine


  !> just force a zero-order hold for the ghost cells, should probably only be used for debugging purposes
  subroutine forceZOH(param)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: param
    integer :: lx1,lx2,lx3

    lx1=size(param,1)-4; lx2=size(param,2)-4; lx3=size(param,3)-4;

!    param(0,:,:,:)=param(1,:,:,:)
!    param(-1,:,:,:)=param(1,:,:,:)
!    param(lx1+1,:,:,:)=param(lx1,:,:,:)
!    param(lx1+2,:,:,:)=param(lx1,:,:,:)
!
!    param(:,0,:,:)=param(:,1,:,:)
!    param(:,-1,:,:)=param(:,1,:,:)
!    param(:,lx2+1,:,:)=param(:,lx2,:,:)
!    param(:,lx2+2,:,:)=param(:,lx2,:,:)

    param(:,:,0,:)=param(:,:,1,:)
    param(:,:,-1,:)=param(:,:,1,:)
    param(:,:,lx3+1,:)=param(:,:,lx3,:)
    param(:,:,lx3+2,:)=param(:,:,lx3,:)
  end subroutine forceZOH


  !> force all primary variables to have a zero-order hold extrapolation in ghost cells for testing purposes
  subroutine forceZOH_all(fluidvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call forceZOH(ns)
    call forceZOH(vs1)
    call forceZOH(vs2)
    call forceZOH(vs3)
    call forceZOH(Ts)
  end subroutine forceZOH_all


  !> permute state variables x1,x2,x3 --> x2,x3,x1.  FIXME: this is now deprecated since we also need to 
  !    swap for density variables when using forestclaw.
  subroutine permute_fluidvars(fluidvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    integer lx1,lx2,lx3,leqn,ix1,ix2,ix3,ieqn
    real(wp), dimension(:,:,:,:), allocatable :: fluidvarsT

    ! size and allocate
    lx1=size(fluidvars,1)-4; lx2=size(fluidvars,2)-4; lx3=size(fluidvars,3)-4; leqn=size(fluidvars,4);
    allocate(fluidvarsT(-1:lx2+2,-1:lx3+2,-1:lx1+2,leqn))

    ! place data in permuted array
    do ieqn=1,leqn
      do ix3=-1,lx3+2
        do ix2=-1,lx2+2
          do ix1=-1,lx1+2
            fluidvarsT(ix2,ix3,ix1,ieqn)=fluidvars(ix1+2,ix2+2,ix3+2,ieqn)
          end do
        end do
      end do
    end do

    ! reshape the permuted array to match original, copy into state variable array
    fluidvars=reshape(fluidvarsT,[lx1+4,lx2+4,lx3+4,leqn])

    ! explicitly clear memory
    deallocate(fluidvarsT)
  end subroutine permute_fluidvars


  ! inverse permute state variables x2,x3,x1 --> x1,x2,x3.  FIXME: also deprecated.  
  subroutine ipermute_fluidvars(fluidvars)
    real(wp), dimension(:,:,:,:), pointer, intent(inout) :: fluidvars
    integer lx1,lx2,lx3,leqn,ix1,ix2,ix3,ieqn
    real(wp), dimension(:,:,:,:), allocatable :: fluidvarsT

    ! size and allocate, sizes of input array will be correct but data will be shaped wrongly
    lx1=size(fluidvars,1)-4; lx2=size(fluidvars,2)-4; lx3=size(fluidvars,3)-4; leqn=size(fluidvars,4);
    allocate(fluidvarsT(-1:lx2+2,-1:lx3+2,-1:lx1+2,leqn))
    fluidvarsT=reshape(fluidvars,[lx2+4,lx3+4,lx1+4,leqn])

    ! place data in permuted array
    do ieqn=1,leqn
      do ix3=-1,lx3+2
        do ix2=-1,lx2+2
          do ix1=-1,lx1+2
            fluidvars(ix1+2,ix2+2,ix3+2,ieqn)=fluidvarsT(ix2,ix3,ix1,ieqn)
          end do
        end do
      end do
    end do

    ! explicitly clear memory
    deallocate(fluidvarsT)
  end subroutine ipermute_fluidvars


  !> Tag for refinement based on location
  subroutine tag4refine_location(x,flagrefine)
    class(curvmesh), intent(in) :: x
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3
    real(wp) :: mlat,mlon,alt

    flagrefine=.false.
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
!          if ( mlon > 209.5 .and. mlon < 210.5 .and. mlat > 28.0 .and. mlat < 29.0 .and. &
!               alt > 80e3 .and. alt < 300e3) then
!            flagrefine=.true.
!            exit
!          end if
          if ( sqrt( ((mlon-210)/3)**2 + (mlat-28.5)**2) < 0.35 .and. &
               alt > 80e3 .and. alt < 300e3) then
            flagrefine=.true.
            exit
          end if
        end do
      end do
    end do
  end subroutine tag4refine_location


  !> Tag for refinement based on perpendicular velocity
  subroutine tag4refine_vperp(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: vperp,vpar

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)


    flagrefine=.false.
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          !vperp=sqrt(vs2(ix1,ix2,ix3,1)**2 + vs3(ix1,ix2,ix3,1)**2)    ! use the major ion species
          !vpar=abs(vs1(ix1,ix2,ix3,1))
          vpar=abs(intvars%atmos%vn1(ix1,ix2,ix3))
          if (alt > 90e3 .and. alt < 350e3 .and. vpar > 25._wp) then     ! more than 50 m/s probably  means something is happening
            flagrefine=.true.
            exit
          end if
        end do
      end do
    end do
  end subroutine tag4refine_vperp


  !> Tag for refinement based on differences within some range
!  subroutine tag4refine_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
!    type(gemini_work) :: intvars
!    logical, intent(inout) :: flagrefine
!    integer :: ix1,ix2,ix3    
!    real(wp) :: mlat,mlon,alt
!    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
!    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
!    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
!    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: minv,maxv,deltav
!    
!    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
!    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!
!    flagrefine=.false.
!    deltav=0.0
!    minv=0.0
!    maxv=0.0
!    do ix3=1,x%lx3
!      do ix2=1,x%lx2
!        do ix1=1,x%lx1
!          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
!          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
!          alt=x%alt(ix1,ix2,ix3)
!          if (sqrt( ((mlon-210)/3)**2 + (mlat-28.5)**2) < 0.35 .and. &
!               alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
!            if (vs1(ix1,ix2,ix3,7) < minv) minv=vs1(ix1,ix2,ix3,7)
!            if (vs1(ix1,ix2,ix3,7) > maxv) maxv=vs1(ix1,ix2,ix3,7)
!          end if
!        end do
!      end do
!    end do
!    deltav=maxv-minv
!
!    if (deltav > 10.0) flagrefine=.true.
!  end subroutine tag4refine_diff

  subroutine tag4refine_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3    
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: minv,maxv,deltav
    
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    flagrefine=.false.
    deltav=0.0
    minv=0.0
    maxv=0.0
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          if (alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
            if (intvars%atmos%vn1(ix1,ix2,ix3) < minv) minv=intvars%atmos%vn1(ix1,ix2,ix3)
            if (intvars%atmos%vn1(ix1,ix2,ix3) > maxv) maxv=intvars%atmos%vn1(ix1,ix2,ix3)
          end if
        end do
      end do
    end do
    deltav=maxv-minv

    if (deltav > 10.0) flagrefine=.true.
  end subroutine tag4refine_diff


  !> Tag for refinement based on differences within some range
  subroutine tag4refine_grad(x,fluidvars,fluidauxvars,electrovars,intvars,flagrefine)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagrefine
    integer :: ix1,ix2,ix3    
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: vtest1,vtest2
    real(wp), dimension(x%lx1,x%lx2,x%lx3) :: gradvn12,gradvn13
    
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

!    vtest1=minval(intvars%atmos%vn1)
!    vtest2=maxval(intvars%atmos%vn1)
!    if (vtest2-vtest1<15._wp) then
!      flagrefine=.false.
!    else
!      flagrefine=.true.
!    end if

    gradvn12(:,:,:)=grad3D2(intvars%atmos%vn1(:,:,:),x,1,x%lx1,1,x%lx2,1,x%lx3)
    gradvn13(:,:,:)=grad3D3(intvars%atmos%vn1(:,:,:),x,1,x%lx1,1,x%lx2,1,x%lx3)

    flagrefine=.false.
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          if (alt > 90e3 .and. alt < 350e3 .and. &
                  (abs(gradvn12(ix1,ix2,ix3)) > 1.25e-3 .or. abs(gradvn13(ix1,ix2,ix3)) > 1.25e-3) ) then
            flagrefine=.true.
            exit
          end if
        end do
      end do
    end do
  end subroutine tag4refine_grad


    !> Tag for coarsening based on differences within some range
!  subroutine tag4coarsening_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagcoarsening)
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
!    type(gemini_work) :: intvars
!    logical, intent(inout) :: flagcoarsening
!    integer :: ix1,ix2,ix3    
!    real(wp) :: mlat,mlon,alt
!    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
!    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
!    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
!    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
!    real(wp) :: minv,maxv,deltav
!    
!    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
!    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
!    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)
!
!    flagcoarsening=.false.
!    deltav=0.0
!    minv=0.0
!    maxv=0.0
!    do ix3=1,x%lx3
!      do ix2=1,x%lx2
!        do ix1=1,x%lx1
!          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
!          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
!          alt=x%alt(ix1,ix2,ix3)
!          if (sqrt( ((mlon-210)/3)**2 + (mlat-28.5)**2) < 0.35 .and. &
!               alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
!            if (vs1(ix1,ix2,ix3,7) < minv) minv=vs1(ix1,ix2,ix3,7)
!            if (vs1(ix1,ix2,ix3,7) > maxv) maxv=vs1(ix1,ix2,ix3,7)
!          end if
!        end do
!      end do
!    end do
!    deltav=maxv-minv
!
!    if (deltav < 5.0) flagcoarsening=.true.
!  end subroutine tag4coarsening_diff

  subroutine tag4coarsening_diff(x,fluidvars,fluidauxvars,electrovars,intvars,flagcoarsening)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars,electrovars
    type(gemini_work) :: intvars
    logical, intent(inout) :: flagcoarsening
    integer :: ix1,ix2,ix3    
    real(wp) :: mlat,mlon,alt
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:), pointer :: E1,E2,E3,J1,J2,J3,Phi
    real(wp) :: minv,maxv,deltav
    
    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)
    call electrovar_pointers(electrovars,E1,E2,E3,J1,J2,J3,Phi)

    flagcoarsening=.false.
    deltav=0.0
    minv=0.0
    maxv=0.0
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          mlon=x%phi(ix1,ix2,ix3)*180.0/pi
          mlat=90.0-x%theta(ix1,ix2,ix3)*180.0/pi
          alt=x%alt(ix1,ix2,ix3)
          if (alt > 80e3 .and. alt < 300e3) then      ! only update min/max v if in region of interest
            if (intvars%atmos%vn1(ix1,ix2,ix3) < minv) minv=intvars%atmos%vn1(ix1,ix2,ix3)
            if (intvars%atmos%vn1(ix1,ix2,ix3) > maxv) maxv=intvars%atmos%vn1(ix1,ix2,ix3)
          end if
        end do
      end do
    end do
    deltav=maxv-minv

    if (deltav < 5.0) flagcoarsening=.true.
  end subroutine tag4coarsening_diff


  !> if refinement is being done it may be advantageous to have refine/interpolate done with drift and temperature
  !    it's easiest to just copy swap existing variables.
  subroutine swap_statevars(fluidvars,fluidauxvars)
    real(wp), dimension(:,:,:,:), pointer, intent(in) :: fluidvars,fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:,:), pointer :: rhovs1,rhoes
    real(wp), dimension(:,:,:), pointer :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom
    real(wp), dimension(:,:,:,:), allocatable :: tmpvar

    call fluidvar_pointers(fluidvars,ns,vs1,vs2,vs3,Ts)
    call fluidauxvar_pointers(fluidauxvars,rhovs1,rhoes,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom)

    !allocate(tmpvar(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2,1:size(ns,4)))
    allocate(tmpvar,mold=rhovs1)
    tmpvar(:,:,:,:)=rhovs1(:,:,:,:)
    rhovs1(:,:,:,:)=vs1(:,:,:,:)
    vs1(:,:,:,:)=tmpvar(:,:,:,:)

    tmpvar(:,:,:,:)=rhoes(:,:,:,:)
    rhoes(:,:,:,:)=Ts(:,:,:,:)
    Ts(:,:,:,:)=tmpvar(:,:,:,:)
    deallocate(tmpvar)
  end subroutine swap_statevars


!  subroutine tag4refine_dist(x,flagrefine)
!    class(curvmesh), intent(in) :: x
!    logical, intent(inout) :: flagrefine
!    real(wp) :: phi1,theta1
!    real(wp), dimension(1:x%lx1,1:x%lx2,1:x%lx3) :: xi,yi,zi
!    real(wp), parameter :: rmag=25e3
!    real(wp) :: r,xi0,yi0
!    integer :: ix1,ix2,ix3
!
!    ! retrieve grid ENU coords
!    call geog2geomag(x%glonctr,x%glatctr,phi1,theta1)
!    call ECEFspher2ENU(x%alt(1:x%lx1,1:x%lx2,1:x%lx3),x%theta(1:x%lx1,1:x%lx2,1:x%lx3), &
!                         x%phi(1:x%lx1,1:x%lx2,1:x%lx3), &
!                         theta1,phi1,xi,yi,zi)
!
!    ! find the average x,y on the end slab so we can control refinement
!    ix1=x%lx1-25
!    xi0=sum(xi(ix1,:,:))/(x%lx2*x%lx3)
!    yi0=sum(yi(ix1,:,:))/(x%lx2*x%lx3)
!
!    ! now check for distance from center to determine refinment
!    flagrefine=.false.
!    do ix3=1,lx3
!      do ix2=1,lx2
!          r=sqrt( (xi(ix1,ix2,ix3)-xi0)**2 + (yi(ix1,ix2,ix3)-yi0)**2)
!          if (r<rmag) then
!            flagrefine=.true.
!            exit
!          end if
!      end do
!    end do
!  end subroutine tag4refine_dist


  !> Call function to retrieve locations from a neutraldata3D_fclaw object
  subroutine get_locationsi_in(intvars,flagallpts,zlims,xlims,ylims,zvals,xvals,yvals,datavals)
    type(gemini_work), intent(inout) :: intvars
    logical, intent(in) :: flagallpts
    real(wp), dimension(2), intent(in) :: zlims,xlims,ylims
    real(wp), dimension(:), pointer, intent(inout) :: zvals,xvals,yvals
    real(wp), dimension(:,:), pointer, intent(inout) :: datavals
    class(neutraldata), pointer :: aperptr

    aperptr=>intvars%atmosperturb    ! apparently select case cannot handle a compound statement
    select type (aperptr)
    class is (neutraldata3D_fclaw)
      call intvars%atmosperturb%get_locationsi(flagallpts,zlims,xlims,ylims,zvals,xvals,yvals,datavals)
    class default
      print*, 'WARNING:  attempted to direct feed data (get) to object of wrong type (not neutraldata3D_fclaw)'
      zvals=>null()
      xvals=>null()
      yvals=>null()
      datavals=>null()
    end select
  end subroutine get_locationsi_in


  !> retrieve a pointer that we can directly copy data to
  subroutine get_datainow_ptr_in(intvars,datavals)
    type(gemini_work), intent(inout) :: intvars
    real(wp), dimension(:,:), pointer, intent(inout) :: datavals
    class(neutraldata), pointer :: aperptr

    aperptr=>intvars%atmosperturb    ! apparently select case cannot handle a compound statement
    select type (aperptr)
    class is (neutraldata3D_fclaw)
      datavals=>intvars%atmosperturb%get_datainow_ptr()
    class default
      print*, 'WARNING:  attempted to direct feed data (get) to object of wrong type (not neutraldata3D_fclaw)'
      datavals=>null();
    end select
  end subroutine get_datainow_ptr_in


  !> Notify object that the data have been placed in its buffer and it needs to copy-out
  subroutine set_datainow_in(intvars)
    type(gemini_work), intent(inout) :: intvars
    class(neutraldata), pointer :: aperptr

    aperptr=>intvars%atmosperturb    ! apparently select case cannot handle a compound statement
    select type (aperptr)
    class is (neutraldata3D_fclaw)
      call intvars%atmosperturb%set_datainow()
    class default
      print*, 'WARNING:  attempted to direct feed data (set) to object of wrong type (not neutraldata3D_fclaw)'
    end select
  end subroutine set_datainow_in


  !> Interface to access trilinear interpolation routines in gemini and interpolate MAGIC data
  subroutine interp3_in(x,y,z,q,aux,xi,yi,zi,qi,auxi)
    real(wp), dimension(:), intent(in) :: x
    real(wp), dimension(:), intent(in) :: y
    real(wp), dimension(:), intent(in) :: z
    real(wp), intent(in) :: xi,yi,zi    ! single points based on forestclaw organization
    real(wp), dimension(:), intent(inout) :: qi,auxi
    real(wp), intent(in), dimension(:,:,:,:) :: q,aux
    integer :: mx,my,mz,meqn,maux,ieqn,iaux
    real(wp), dimension(1) :: xiarr,yiarr,ziarr,qiarr,auxiarr

    ! sizes
    mx=size(q,1); my=size(q,2); mz=size(q,3);
    meqn=size(q,4); maux=size(aux,4);
    xiarr=[xi]; yiarr=[yi]; ziarr=[zi];

    ! interpolate
    do ieqn=1,meqn
      qiarr=interp3(x,y,z,q(:,:,:,ieqn),xiarr,yiarr,ziarr)     !note that the MAGIC coordinatees are permuted x,y,z
      qi(ieqn)=qiarr(1)
    end do
    do iaux=1,maux
      auxiarr=interp3(x,y,z,aux(:,:,:,iaux),xiarr,yiarr,ziarr)
      auxi(iaux)=auxiarr(1)
    end do
  end subroutine interp3_in


  !> Interface to access trilinear interpolation routines in gemini and interpolate MAGIC data
  subroutine interp2_in(x,y,q,aux,rhoi,zi,qi,auxi)
    real(wp), dimension(:), intent(in) :: x
    real(wp), dimension(:), intent(in) :: y
    real(wp), intent(in) :: rhoi,zi    ! single points based on forestclaw organization
    real(wp), dimension(:), intent(inout) :: qi,auxi
    real(wp), intent(in), dimension(:,:,:) :: q,aux
    integer :: mx,my,meqn,maux,ieqn,iaux
    real(wp), dimension(1) :: xiarr,yiarr,qiarr,auxiarr

    ! sizes
    !mx=size(q,1); my=size(q,2);
    meqn=size(q,3); maux=size(aux,3);
    xiarr=[rhoi]; yiarr=[zi];

    ! interpolate
    do ieqn=1,meqn
      qiarr=interp2(x,y,q(:,:,ieqn),xiarr,yiarr)     !note that the MAGIC coordinatees are permuted x,y,z
      qi(ieqn)=qiarr(1)
    end do
    do iaux=1,maux
      auxiarr=interp2(x,y,aux(:,:,iaux),xiarr,yiarr)
      auxi(iaux)=auxiarr(1)
    end do
  end subroutine interp2_in


  !> increment date and time arrays, this is superfluous but trying to keep outward facing function calls here.
  subroutine dateinc_in(dt,ymd,UTsec)
    real(wp), intent(in) :: dt
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec

    call dateinc(dt,ymd,UTsec)

    !print*, 'Date updated to:  ',ymd,UTsec
  end subroutine dateinc_in
end module gemini3d
