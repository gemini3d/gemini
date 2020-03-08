include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

# for CMake >= 3.17 with  cmake -G "Ninja Multi-Config"
set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "Build type selections" FORCE)

# Do these before compiler options so options don't goof up finding
# === OpenMP
# optional, for possible MUMPS speedup
if(openmp)
  find_package(OpenMP COMPONENTS C Fortran)
endif()
# === MPI
# MPI is used throughout Gemini
find_package(MPI REQUIRED COMPONENTS Fortran)

# === compiler setup
# feel free to add more compiler_*.cmake
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_intel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_gnu.cmake)
endif()
