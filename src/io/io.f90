module io
!! HANDLES INPUT AND OUTPUT OF PLASMA STATE PARAMETERS (NOT GRID INPUTS)
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use gemini3d_config, only : gemini_cfg
use phys_consts, only : kB,ms,pi,lsp,wp,lwave, comp_lvl
use mpimod, only: mpi_cfg, tag=>gemini_mpi
use grid, only : gridflag,lx1,lx2,lx3,lx2all, lx3all
use gemini_work_def, only: gemini_work

implicit none (type, external)

private
public :: create_outdir, &
  input_plasma, output_plasma, input_plasma_currents, &
  create_outdir_mag, output_magfields, &
  output_aur, output_cond, &
  find_milestone


interface !< aurora.f90
  module subroutine output_aur(outdir,flagglow,ymd,UTsec,iver, out_format)
    character(*), intent(in) :: outdir, out_format
    integer, intent(in) :: flagglow
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:), intent(in) :: iver
  end subroutine output_aur
end interface


interface !< cond.f90
  module subroutine output_cond(outdir, ymd, UTsec, sig0, sigP, sigH, out_format)
    character(*), intent(in) :: outdir, out_format
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:), intent(in) :: sig0, sigP, sigH
  end subroutine output_cond
end interface


interface !< mag.f90
  module subroutine create_outdir_mag(outdir,fieldpointfile)
    character(*), intent(in) :: outdir
    character(*), intent(in) :: fieldpointfile
  end subroutine create_outdir_mag

  module subroutine output_magfields(outdir,ymd,UTsec,Br,Btheta,Bphi, out_format)
    character(*), intent(in) :: outdir, out_format
    integer, intent(in) :: ymd(3)
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:), intent(in)  :: Br,Btheta,Bphi
  end subroutine output_magfields
end interface


interface !< plasma.f90
  module subroutine input_plasma(out_dir, x1,x2,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
    character(*), intent(in) :: out_dir
    real(wp), dimension(-1:), intent(in) :: x1, x2, x3all
    character(*), intent(in) :: indatsize, indatfile
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,vs1,Ts
    !! intent(out)
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phi
    !! intent(out)
    real(wp), dimension(:,:,:), pointer, intent(inout) :: Phiall
    !! intent(out)
  end subroutine input_plasma

  module subroutine input_plasma_currents(outdir,out_format,flagoutput,ymd,UTsec,J1,J2,J3)
    character(*), intent(in) :: outdir, out_format
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J1,J2,J3
    !! intent(out)
  end subroutine input_plasma_currents

  module subroutine output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3,out_format,intvars)
    character(*), intent(in) :: outdir, out_format
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
    real(wp), dimension(:,:,:), pointer, intent(inout) :: Phiall     !these jokers may not be allocated, but this is allowed as of f2003
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: J1,J2,J3
    type(gemini_work), intent(in) :: intvars
  end subroutine output_plasma
end interface


interface !< output.f90
  module subroutine create_outdir(cfg)
    class(gemini_cfg), intent(in) :: cfg
  end subroutine create_outdir
end interface


interface !< milestone.f90
  module subroutine find_milestone(cfg, tmile,ymdmile,UTsecmile,filemile)
    class(gemini_cfg), intent(in) :: cfg
    real(wp), intent(out) :: tmile
    !! time elapsed from beginning of the simulation to the milestone
    integer, dimension(3), intent(out) :: ymdmile
    real(wp), intent(out) :: UTsecmile
    character(:), allocatable, intent(out) :: filemile
  end subroutine find_milestone
end interface

end module io
