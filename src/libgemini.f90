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
!!   that does not involve mpi or mpi-dependent modules.  
module gemini3d

use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float
use gemini_cli, only : cli
use gemini_init, only : find_config, check_input_files
use phys_consts, only: wp,debug,lnchem,lwave,lsp
use meshobj, only: curvmesh
use config, only: gemini_cfg
use collisions, only: conductivities
use pathlib, only : expanduser
use grid, only: grid_size,x,lx1,lx2,lx3
use config, only : read_configfile,cfg

implicit none (type, external)
private
public :: c_params, cli_config_gridsize, gemini_alloc, gemini_dealloc, init_avgparms_Bfield

!> type for passing C-like parameters between program units
type, bind(C) :: c_params
  !! this MUST match gemini3d.h and libgemini.f90 exactly including order
  logical(c_bool) :: fortran_cli, debug, dryrun
  character(kind=c_char) :: out_dir(1000)
  !! .ini [base]
  integer(c_int) :: ymd(3)
  real(kind=c_float) :: UTsec0, tdur, dtout, activ(3), tcfl, Teinf
  !! .ini
end type c_params

!> libgem_utils.f90
interface
  module subroutine cli_config_gridsize(p,lid2in,lid3in)
    type(c_params), intent(in) :: p
    integer, intent(inout) :: lid2in,lid3in
  end subroutine cli_config_gridsize
  module subroutine gemini_alloc(ns,vs1,vs2,vs3,Ts,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom, &
                                    E1,E2,E3,J1,J2,J3,Phi,nn,Tn,vn1,vn2,vn3,iver)
    real(wp), dimension(:,:,:,:), allocatable, intent(inout) :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom,E1,E2,E3,J1,J2,J3,Phi
    real(wp), dimension(:,:,:,:), allocatable, intent(inout) :: nn
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: Tn,vn1,vn2,vn3
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: iver
  end subroutine gemini_alloc
  module subroutine gemini_dealloc(ns,vs1,vs2,vs3,Ts,rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom,& 
                                      E1,E2,E3,J1,J2,J3,Phi,nn,Tn,vn1,vn2,vn3,iver)
    real(wp), dimension(:,:,:,:), allocatable, intent(inout) :: ns,vs1,vs2,vs3,Ts
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: rhov2,rhov3,B1,B2,B3,v1,v2,v3,rhom,E1,E2,E3,J1,J2,J3,Phi
    real(wp), dimension(:,:,:,:), allocatable, intent(inout) :: nn
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: Tn,vn1,vn2,vn3
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: iver
  end subroutine  gemini_dealloc
  module subroutine init_avgparms_Bfield(rhov2,rhov3,v2,v3,B2,B3,B1)
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: rhov2,rhov3,v2,v3,B1,B2,B3
  end subroutine init_avgparms_Bfield
end interface

contains

!> NOTHING other than submodules :)

end module gemini3d
