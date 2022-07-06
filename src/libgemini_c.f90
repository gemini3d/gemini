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

!! This module contains C/CXX wrappers for functions in libgemini.
!! These routines match those in libgemini.f90 and are
!! principally meant to convert the C pointers to various data objects into fortran pointers.
!! The grid is a class pointer (pointer to polymorphic object).
!! Other polymorphic objects (neutraldata, etc.) are kept in a static
!! derived type (intvars::gemini_work) and don't need to be passes as class pointers.

module gemini3d_C

use, intrinsic :: iso_c_binding, only : c_int, c_bool, c_loc, c_null_ptr, c_ptr, c_f_pointer, wp => C_DOUBLE
use phys_consts, only: lnchem,lwave,lsp
use grid, only: lx1,lx2,lx3, detect_gridtype
use meshobj, only: curvmesh
use meshobj_cart, only: cartmesh
use meshobj_dipole, only: dipolemesh
use precipdataobj, only: precipdata
use efielddataobj, only: efielddata
use neutraldataobj, only: neutraldata
use gemini3d_config, only: gemini_cfg
use gemini3d, only: c_params, cli_config_gridsize, init_precipinput_in, msisinit_in, &
            set_start_values, init_neutralBG_in, set_update_cadence, neutral_atmos_winds, get_solar_indices, &
            v12rhov1_in, T2rhoe_in, interface_vels_allspec_in, sweep3_allparams_in, &
            sweep1_allparams_in, sweep2_allparams_in, &
            rhov12v1_in, VNRicht_artvisc_in, compression_in, rhoe2T_in, clean_param_in, &
            energy_diffusion_in, source_loss_allparams_in, &
            dateinc_in, get_subgrid_size,get_fullgrid_size,get_config_vars, get_species_size, fluidvar_pointers, &
            fluidauxvar_pointers, electrovar_pointers, gemini_work, gemini_alloc_nodouble, gemini_dealloc_nodouble, &
            interp_file2subgrid_in,grid_from_extents_in,read_fullsize_gridcenter_in, setv2v3

implicit none (type, external)

public

contains
  !> set fortran object pointer dynamic type to what is indicated in objtype.  Convert C pointer using
  !>    declared static types (c_f_pointer will not work on a polymorphic object).
  function set_gridpointer_dyntype(xtype,xC) result(x)
    type(c_ptr), intent(in) :: xC
    integer(C_INT), intent(in) :: xtype
    class(curvmesh), pointer :: x
    type(cartmesh), pointer :: xcart
    type(dipolemesh), pointer :: xdipole

    select case (xtype)
      case (1)
        call c_f_pointer(xC,xcart)
        x=>xcart
      case (2)
        call c_f_pointer(xC,xdipole)
        x=>xdipole
      case default
        error stop 'unable to identify object type during conversion from C to fortran class pointer'
    end select
  end function set_gridpointer_dyntype


  !> basic command line and grid size determination
  subroutine cli_config_gridsize_C(p,lid2in,lid3in,cfgC) bind(C, name='cli_config_gridsize_C')
    type(c_params), intent(in) :: p
    integer(C_INT), intent(inout) :: lid2in,lid3in
    type(c_ptr), intent(inout) :: cfgC
    type(gemini_cfg), pointer :: cfg

    allocate(cfg)
    call cli_config_gridsize(p,lid2in,lid3in,cfg)
    cfgC = c_loc(cfg)
  end subroutine cli_config_gridsize_C


  !> return some data from cfg that is needed in the main program
  subroutine get_config_vars_C(cfgC,flagneuBG,flagdneu,dtneuBG,dtneu) bind(C, name='get_config_vars_C')
    type(c_ptr), intent(in) :: cfgC
    logical(C_BOOL), intent(inout) :: flagneuBG
    integer(C_INT), intent(inout) :: flagdneu
    real(wp), intent(inout) :: dtneuBG,dtneu

    type(gemini_cfg), pointer :: cfg
    logical :: neuBG

    neuBG = flagneuBG

    call c_f_pointer(cfgC,cfg)
    call get_config_vars(cfg, neuBG, flagdneu,dtneuBG,dtneu)

    flagneuBG = neuBG
  end subroutine get_config_vars_C


  !> returns the subgrid sizes *** stored in the grid module ***
  subroutine get_subgrid_size_C(lx1out,lx2out,lx3out) bind(C, name='get_subgrid_size_C')
    integer(C_INT), intent(inout) :: lx1out,lx2out,lx3out

    call get_subgrid_size(lx1out,lx2out,lx3out)
  end subroutine get_subgrid_size_C


  !> return full grid extents *** stored in the grid module ***
  subroutine get_fullgrid_size_C(lx1out,lx2allout,lx3allout) bind(C, name='get_fullgrid_size_C')
    integer(C_INT), intent(inout) :: lx1out,lx2allout,lx3allout

    call get_fullgrid_size(lx1out, lx2allout, lx3allout)
  end subroutine get_fullgrid_size_C


  !> return number of species *** from phys_consts module ***
  subroutine get_species_size_C(lspout) bind(C, name='get_species_size_C')
    integer(C_INT), intent(inout) :: lspout

    call get_species_size(lspout)
  end subroutine get_species_size_C


  !> force a value for the grid drifts
  subroutine setv2v3_C(v2gridin,v3gridin) bind(C, name='setv2v3_C')
    real(wp), intent(in) :: v2gridin,v3gridin

    call setv2v3(v2gridin,v3gridin)
  end subroutine setv2v3_C


  !> allocate space for gemini state variables, bind pointers to blocks of memory
  subroutine gemini_alloc_C(cfgC,intvarsC) bind(C, name='gemini_alloc_C')
    type(c_ptr), intent(in) :: cfgC
    type(c_ptr), intent(inout) :: intvarsC
    type(gemini_cfg), pointer :: cfg
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    allocate(intvars)
    call gemini_alloc_nodouble(cfg,intvars)
    intvarsC=c_loc(intvars)
  end subroutine gemini_alloc_C


  !> deallocate state variables
  subroutine gemini_dealloc_C(cfgC,intvarsC) bind(C, name='gemini_dealloc_C')
    type(c_ptr), intent(in) :: cfgC
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    call c_f_pointer(intvarsC,intvars)

    !> there are issues with allocating primitives variables (doubles/ints) and then deallocating
    !    when passed back and forth with C so only deallocate the derived types
    call gemini_dealloc_nodouble(cfg,intvars)
  end subroutine gemini_dealloc_C


  !> C wrapper for procedure to get the center location of the grid from its input file
  subroutine read_fullsize_gridcenter_C(cfgC) bind(C,name='set_fullsize_gridcenter_C')
    type(c_ptr), intent(in) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call read_fullsize_gridcenter_in(cfg)
  end subroutine read_fullsize_gridcenter_C


  !> C wrapper for procedure to compute a grid object given extents and fullgrid reference point
  subroutine grid_from_extents_C(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,xtype,xC) bind(C,name='grid_from_extents_C')
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    integer, intent(inout) :: xtype
    type(c_ptr), intent(inout) :: xC
    class(curvmesh), pointer :: x

    x=>set_gridpointer_dyntype(xtype,xC)
    call grid_from_extents_in(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)  
    ! as an extra step we need to also assign a type to the grid
    xtype=detect_gridtype(x%x1,x%x2,x%x3)
  end subroutine grid_from_extents_C


  !> C wrapper for procedure that reads full data from input file and interpolates to loca worker subgrid
  subroutine interp_file2subgrid_C(cfgC,xtype,xC,fluidvarsC,electrovarsC) bind(C,name='interp_file2subgrid_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: electrovarsC
    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call interp_file2subgrid_in(cfg,x,fluidvars,electrovars)
  end subroutine interp_file2subgrid_C


  !> set start values for some variables.
  !    some case is required here because the state variable pointers are mapped;
  !    however, note that the lbound and ubound have not been set since arrays
  !    are not passed through as dummy args
  !    with specific ubound so that we need to use intrinsic calls to make sure we fill
  !    computational cells (not ghost)
  subroutine set_start_values_C(it,t,tout,tglowout,tneuBG,xtype,xC,fluidauxvarsC) bind(C, name='set_start_values_C')
    integer(C_INT), intent(inout) :: it
    real(wp), intent(inout) :: t,tout,tglowout,tneuBG
    type(c_ptr), intent(inout) :: xC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    integer(C_INT), intent(in) :: xtype
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call set_start_values(it,t,tout,tglowout,tneuBG,x,fluidauxvars)
  end subroutine set_start_values_C


  !> Wrapper for initialization of electron precipitation data
  subroutine init_precipinput_C(cfgC,xtype,xC,dt,t,ymd,UTsec,intvarsC) bind(C, name='init_precipinput_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(intvarsC,intvars)
    call init_precipinput_in(cfg,x,dt,t,ymd,UTsec,intvars)
  end subroutine init_precipinput_C


  !> initialization procedure needed for MSIS 2.0
  subroutine msisinit_C(cfgC) bind(C, name='msisinit_C')
    type(c_ptr), intent(in) :: cfgC
    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call msisinit_in(cfg)
  end subroutine msisinit_C


  !> call to initialize the neutral background information
  subroutine init_neutralBG_C(cfgC,xtype,xC,dt,t,ymd,UTsec,intvarsC) bind(C, name='init_neutralBG_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x    ! so neutral module can deallocate unit vectors once used...
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype,xC)
    call c_f_pointer(intvarsC,intvars)
    call init_neutralBG_in(cfg,x,dt,t,ymd,UTsec,intvars)
  end subroutine init_neutralBG_C


  !> set update cadence for printing out diagnostic information during simulation
  subroutine set_update_cadence_C(iupdate) bind(C, name='set_update_cadence_C')
    integer(C_INT), intent(inout) :: iupdate

    call set_update_cadence(iupdate)
  end subroutine set_update_cadence_C


  !> compute background neutral density, temperature, and wind
  subroutine neutral_atmos_winds_C(cfgC,xtype,xC,ymd,UTsec,intvarsC) bind(C, name='neutral_atmos_winds_C')
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    type(c_ptr), intent(inout) :: intvarsC

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC,cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(intvarsC,intvars)
    call neutral_atmos_winds(cfg,x,ymd,UTsec,intvars)
  end subroutine neutral_atmos_winds_C


  !> get solar indices from cfg struct
  subroutine get_solar_indices_C(cfgC,f107,f107a) bind(C, name='get_solar_indices_C')
    type(c_ptr), intent(in) :: cfgC
    real(wp), intent(inout) :: f107,f107a

    type(gemini_cfg), pointer :: cfg

    call c_f_pointer(cfgC,cfg)
    call get_solar_indices(cfg,f107,f107a)
  end subroutine get_solar_indices_C


  !> convert velocity to momentum density
  subroutine v12rhov1_C(fluidvarsC,fluidauxvarsC) bind(C,name='v12rhov1_C')
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call v12rhov1_in(fluidvars,fluidauxvars)
  end subroutine v12rhov1_C


  !> convert temperature to specific internal energy density
  subroutine T2rhoe_C(fluidvarsC,fluidauxvarsC) bind(C, name='T2rhoe_C')
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call T2rhoe_in(fluidvars,fluidauxvars)
  end subroutine T2rhoe_C


  !> compute interface velocities once haloing has been done
  subroutine interface_vels_allspec_C(fluidvarsC,intvarsC,lsp) bind(C, name='interface_vels_allspec_C')
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: lsp

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call interface_vels_allspec_in(fluidvars,intvars,lsp)
  end subroutine interface_vels_allspec_C


  !> functions for sweeping advection
  subroutine sweep3_allparams_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name='sweep3_allparams_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC, intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep3_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep3_allparams_C

  subroutine sweep1_allparams_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name='sweep1_allparams_C')
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep1_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep1_allparams_C

  subroutine sweep2_allparams_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="sweep2_allparams_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call sweep2_allparams_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine sweep2_allparams_C


  !> conversion of momentum density to velocity
  subroutine rhov12v1_C(fluidvarsC, fluidauxvarsC) bind(C, name="rhov12v1_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call rhov12v1_in(fluidvars,fluidauxvars)
  end subroutine rhov12v1_C


  !> compute artifical viscosity
  subroutine VNRicht_artvisc_C(fluidvarsC,intvarsC) bind(C, name="VNRicht_artvisc_C")
    type(c_ptr), intent(in) :: fluidvarsC
    type(c_ptr), intent(inout) :: intvarsC
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(intvarsC,intvars)
    call VNRicht_artvisc_in(fluidvars,intvars)
  end subroutine VNRicht_artvisc_C


  !> compression substep for fluid solve
  subroutine compression_C(fluidvarsC,fluidauxvarsC,intvarsC,xtype,xC,dt) bind(C, name="compression_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    real(wp), intent(in) :: dt

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call c_f_pointer(intvarsC,intvars)
    x=>set_gridpointer_dyntype(xtype, xC)
    call compression_in(fluidvars,fluidauxvars,intvars,x,dt)
  end subroutine compression_C


  !> convert specific internal energy density into temperature
  subroutine rhoe2T_C(fluidvarsC,fluidauxvarsC) bind(C, name="rhoe2T_C")
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: fluidauxvarsC

    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars

    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp+9)])
    call rhoe2T_in(fluidvars,fluidauxvars)
  end subroutine rhoe2T_C


  !> deal with null cell solutions
  subroutine clean_param_C(iparm,xtype,xC,fluidvarsC) bind(C, name="clean_param_C")
    integer(C_INT), intent(in) :: iparm
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(in) :: fluidvarsC

    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars

    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call clean_param_in(iparm,x,fluidvars)
  end subroutine clean_param_C


  !> diffusion of energy
  subroutine energy_diffusion_C(cfgC,xtype,xC,fluidvarsC,electrovarsC,intvarsC,dt) bind(C, name="energy_diffusion_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt

    type(gemini_cfg), pointer :: cfg
    class(curvmesh), pointer :: x
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)
    call energy_diffusion_in(cfg,x,fluidvars,electrovars,intvars,dt)
  end subroutine energy_diffusion_C


  !> source/loss numerical solutions
  subroutine source_loss_allparams_C(cfgC,fluidvarsC,fluidauxvarsC,electrovarsC,intvarsC,xtype,xC,dt,t,ymd, &
                                        UTsec,f107a,f107,first,gavg,Tninf) bind(C, name="source_loss_allparams_C")
    type(c_ptr), intent(in) :: cfgC
    integer(C_INT), intent(in) :: xtype
    type(c_ptr), intent(in) :: xC
    type(c_ptr), intent(inout) :: fluidvarsC
    type(c_ptr), intent(inout) :: fluidauxvarsC
    type(c_ptr), intent(in) :: electrovarsC
    type(c_ptr), intent(in) :: intvarsC
    real(wp), intent(in) :: dt,t
    integer(C_INT), dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), intent(in) :: f107a,f107
    logical(C_BOOL), intent(in) :: first
    real(wp), intent(in) :: gavg,Tninf

    type(gemini_cfg), pointer :: cfg
    real(wp), dimension(:,:,:,:), pointer :: fluidvars
    real(wp), dimension(:,:,:,:), pointer :: fluidauxvars
    real(wp), dimension(:,:,:,:), pointer :: electrovars
    type(gemini_work), pointer :: intvars
    class(curvmesh), pointer :: x

    call c_f_pointer(cfgC, cfg)
    x=>set_gridpointer_dyntype(xtype, xC)
    call c_f_pointer(fluidvarsC,fluidvars,[(lx1+4),(lx2+4),(lx3+4),(5*lsp)])
    call c_f_pointer(fluidauxvarsC,fluidauxvars,[(lx1+4),(lx2+4),(lx3+4),(2*lsp)+9])
    call c_f_pointer(electrovarsC,electrovars,[(lx1+4),(lx2+4),(lx3+4),7])
    call c_f_pointer(intvarsC,intvars)
    call source_loss_allparams_in(cfg,fluidvars,fluidauxvars,electrovars,intvars,x,dt,t,ymd, &
                                        UTsec,f107a,f107,logical(first),gavg,Tninf)
  end subroutine source_loss_allparams_C


  !> increment date and time arrays, this is superfluous but trying to keep outward facing function calls here.
  subroutine dateinc_C(dt,ymd,UTsec) bind(C, name="dateinc_C")
    real(wp), intent(in) :: dt
    integer(C_INT), dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec

    call dateinc_in(dt,ymd,UTsec)
  end subroutine dateinc_C
end module gemini3d_C
