set_property(DIRECTORY PROPERTY  LABELS potential)

# make mumps usage precision-agnostic
# this is needed because Mumps uses Fortran 95 syntax instead of Fortran 2003 syntax
configure_file(mumps_ifc.in.f90 mumps_ifc.f90 @ONLY)
add_library(mumps_interface OBJECT ${CMAKE_CURRENT_BINARY_DIR}/mumps_ifc.f90)
target_link_libraries(mumps_interface PRIVATE MUMPS::MUMPS)

add_library(PDEelliptic OBJECT
PDEelliptic.f90 elliptic2d.f90 elliptic3d.f90)
target_link_libraries(PDEelliptic PRIVATE mpimod const mumps_interface MUMPS::MUMPS)

add_library(potential OBJECT
potential_comm_mumps.f90 potential_root.f90 potential_worker.f90
potential_mumps.f90 potential2d.f90
boundary_conditions/potentialBCs_mumps.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(potential PRIVATE PDEelliptic calculus collision meshobj efielddata interp timeutils mpimod reader const
inputdata efielddata grid gemini_work_def 
MPI::MPI_Fortran
)

add_library(potential_nompi OBJECT potential_nompi.f90)
target_link_libraries(potential_nompi PRIVATE calculus meshobj const grid)

