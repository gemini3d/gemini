module gemini_work_def

use phys_consts, only: wp
use precipdataobj, only: precipdata
use efielddataobj, only: efielddata
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use neutraldata3Dobj_fclaw, only: neutraldata3D_fclaw
use neutral, only: neutral_info

!> type encapsulating internal arrays and parameters needed by gemini.  This is basically a catch-all for any data
!    in a gemini instance that is needed to advance the solution that must be passed into numerical procedures BUt
!    doesn't conform to simple array shapes.
type gemini_work
  real(wp), dimension(:,:,:), pointer :: Phiall=>null()    !! full-grid potential solution.  To store previous time step value
  real(wp), dimension(:,:,:), pointer :: iver    !! integrated volume emission rate of aurora calculated by GLOW

  !> Other variables used by the fluid solvers
  real(wp), dimension(:,:,:,:), pointer :: vs1i
  real(wp), dimension(:,:,:,:), pointer :: vs2i
  real(wp), dimension(:,:,:,:), pointer :: vs3i
  real(wp), dimension(:,:,:,:), pointer :: Q    ! artificial viscosity

  !> Neutral information for top-level gemini program
  type(neutral_info), pointer :: atmos

  !> Inputdata objects that are needed for each subgrid
  type(precipdata), pointer :: eprecip=>null()
  type(efielddata), pointer :: efield=>null()
  class(neutraldata), pointer :: atmosperturb=>null()   ! not associated by default and may never be associated

  !> User can add any other parameters they want to pass around into this type
  real(wp), dimension(:,:,:), pointer :: sigP=>null()
  real(wp), dimension(:,:,:), pointer :: sigH=>null()
end type gemini_work

contains

end module gemini_work_def
