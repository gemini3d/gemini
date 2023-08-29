submodule (io) plasma_output

use mpimod, only : gather_recv, gather_send

implicit none (type, external)

interface ! plasma_output_*.f90
  module subroutine output_root_stream_mpi_hdf5(outdir,flagoutput,ymd,UTsec,v2avgall,v3avgall,nsall,vs1all,Tsall, &
                                                Phiall,J1all,J2all,J3all,neall,v1avgall,Tavgall,Teall,sigPall,sigHall)
    character(*), intent(in) :: outdir
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:), intent(in) :: v2avgall,v3avgall
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: nsall,vs1all,Tsall
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: Phiall   ! okay to have ghost cells b/c already resides on root.
    real(wp), dimension(1:,1:,1:), intent(in) :: J1all,J2all,J3all   ! tricky/confusing - J1,2,3 have ghost cells but these do not!
    real(wp), dimension(:,:,:), intent(in) :: neall,v1avgall,Tavgall,Teall
    real(wp), dimension(:,:,:), intent(in) :: sigPall,sigHall
  end subroutine
end interface

contains
  module procedure output_plasma
    ! subroutine output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
    !! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
    !! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
    !! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)
    
    if (mpi_cfg%myid == 0) then
      call output_root_stream_mpi(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3,out_format,intvars)
    else
      call output_workers_mpi(vs2,vs3,ns,vs1,Ts,J1,J2,J3,intvars)
    end if
  end procedure output_plasma


  subroutine output_workers_mpi(vs2,vs3,ns,vs1,Ts,J1,J2,J3,intvars)
    !------------------------------------------------------------
    !-------SEND COMPLETE DATA FROM WORKERS TO ROOT PROCESS FOR OUTPUT.
    !-------STATE VARS ARE EXPECTED TO INCLUDE GHOST CELLS
    !------- This is the same regardless of what type of output is
    !------- being done.
    !------------------------------------------------------------
        real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: J1,J2,J3
    type(gemini_work), intent(in) :: intvars
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: v2avg,v3avg
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: Jtmp     ! contiguous buffer for MPI
    
    
    !ONLY AVERAGE DRIFTS PERP TO B NEEDED FOR OUTPUT
    v2avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs2(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v2avg=v2avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
    v3avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs3(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v3avg=v3avg/ns(1:lx1,1:lx2,1:lx3,lsp)
    
    
    !SEND MY GRID DATA TO THE ROOT PROCESS
    call gather_send(v2avg,tag%v2)
    call gather_send(v3avg,tag%v3)
    call gather_send(ns,tag%ns)
    call gather_send(vs1,tag%vs1)
    call gather_send(Ts,tag%Ts)
    
    
    !------- SEND ELECTRODYNAMIC PARAMETERS TO ROOT
    Jtmp=J1(1:lx1,1:lx2,1:lx3)
    call gather_send(Jtmp,tag%J1)
    Jtmp=J2(1:lx1,1:lx2,1:lx3)
    call gather_send(Jtmp,tag%J2)
    Jtmp=J3(1:lx1,1:lx2,1:lx3)
    call gather_send(Jtmp,tag%J3)

    ! the user can add gather operations for their various custom variables stored in intvars
    Jtmp=intvars%sigP(1:lx1,1:lx2,1:lx3)
    call gather_send(Jtmp,tag%sigP)
    Jtmp=intvars%sigH(1:lx1,1:lx2,1:lx3)
    call gather_send(Jtmp,tag%sigH)
  end subroutine output_workers_mpi
  
  
  subroutine output_root_stream_mpi(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3,out_format,intvars)
    !------------------------------------------------------------
    !------- Root needs to gather data and pass to subroutine to
    !------- write to disk in the appropriate format.
    !------------------------------------------------------------
    character(*), intent(in) :: outdir
    character(*), intent(in) :: out_format
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: Phiall
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: J1,J2,J3
    type(gemini_work), intent(in) :: intvars
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: v2avg,v3avg
    real(wp), dimension(-1:lx1+2,-1:lx2all+2,-1:lx3all+2,1:lsp) :: nsall,vs1all,Tsall
    real(wp), dimension(1:lx1,1:lx2all,1:lx3all) :: v2avgall,v3avgall,v1avgall,Tavgall,neall,Teall
    real(wp), dimension(1:lx1,1:lx2,1:lx3) :: Jtmp    ! continguous buffer for variables to be message-passed
    real(wp), dimension(1:lx1,1:lx2all,1:lx3all) :: J1all,J2all,J3all,sigPall,sigHall
    
    
    print *, 'System sizes according to Phiall:  ',lx1,lx2all,lx3all
    !ONLY AVERAGE DRIFTS PERP TO B NEEDED FOR OUTPUT
    v2avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs2(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v2avg=v2avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
    v3avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs3(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v3avg=v3avg/ns(1:lx1,1:lx2,1:lx3,lsp)
    
    !GET THE SUBGRID DATA FORM THE WORKERS
    call gather_recv(v2avg,tag%v2,v2avgall)
    call gather_recv(v3avg,tag%v3,v3avgall)
    call gather_recv(ns,tag%ns,nsall)
    call gather_recv(vs1,tag%vs1,vs1all)
    call gather_recv(Ts,tag%Ts,Tsall)
    
    !> RADD--- NEED TO ALSO GATHER FULL GRID ELECTRODYANMICS PARAMETERS FROM WORKERS
    Jtmp=J1(1:lx1,1:lx2,1:lx3)
    call gather_recv(Jtmp,tag%J1,J1all)
    Jtmp=J2(1:lx1,1:lx2,1:lx3)
    call gather_recv(Jtmp,tag%J2,J2all)
    Jtmp=J3(1:lx1,1:lx2,1:lx3)
    call gather_recv(Jtmp,tag%J3,J3all)
    
    !COMPUTE AVERAGE VALUES FOR ION PLASMA PARAMETERS
    !> possible bottleneck; should have workers help?
    !> also only compute these if they are actually being output
    if (flagoutput==2 .or. flagoutput==3) then
      neall=nsall(1:lx1,1:lx2all,1:lx3all,lsp)
    end if
    if (flagoutput==2) then
      v1avgall=sum(nsall(1:lx1,1:lx2all,1:lx3all,1:lsp-1)*vs1all(1:lx1,1:lx2all,1:lx3all,1:lsp-1),4)
      v1avgall=v1avgall/nsall(1:lx1,1:lx2all,1:lx3all,lsp)    !compute averages for output.
      Tavgall=sum(nsall(1:lx1,1:lx2all,1:lx3all,1:lsp-1)*Tsall(1:lx1,1:lx2all,1:lx3all,1:lsp-1),4)
      Tavgall=Tavgall/nsall(1:lx1,1:lx2all,1:lx3all,lsp)    !compute averages for output.
      Teall=Tsall(1:lx1,1:lx2all,1:lx3all,lsp)
    end if

    ! here the user can add any output custom variables stored in intvars, call signature of output below needs to be adapted to
    !   include all arrays user wants in output file
    Jtmp=intvars%sigP(1:lx1,1:lx2,1:lx3)
    call gather_recv(Jtmp,tag%sigP,sigPall)
    Jtmp=intvars%sigH(1:lx1,1:lx2,1:lx3)
    call gather_recv(Jtmp,tag%sigH,sigHall)
    
    !> Now figure out which type of file we write to
    select case (out_format)
      case ('h5')
        call output_root_stream_mpi_hdf5(outdir,flagoutput,ymd,UTsec,v2avgall,v3avgall,nsall,vs1all,Tsall, &
                                           Phiall,J1all,J2all,J3all,neall,v1avgall,Tavgall,Teall,sigPall,sigHall)
      case default
        error stop 'plasma_output:output_root_stream_api: unknown format' // out_format
    end select
  end subroutine output_root_stream_mpi
end submodule plasma_output
