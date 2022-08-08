# NOTE: don't use -march=native as GCC doesn't support all CPU arches with that option.
# add_compile_options(-mtune=native)

add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fimplicit-none>)

# --- IMPORTANT
add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-Werror=array-bounds;-fcheck=all>")
# --- IMPORTANT: options help trap array indexing/bounds errors at runtime


if(dev)
  add_compile_options(-Wall -Wextra)
  # -Wpedantic makes too many false positives
else(dev)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument;-Wno-unused-variable;-Wno-unused-function>")
endif(dev)

# avoid backtrace that's unusable without -g
add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-fno-backtrace>)

# Wdo-subscript is known to warn on obvious non-problems
check_fortran_compiler_flag(-Wdo-subscript dosubflag)
if(dosubflag)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Wno-do-subscript>)
endif()



# add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,RelWithDebInfo>>:-ffpe-trap=invalid,zero,overflow>")#,underflow)

# lot of spurious warnings on allocatable scalar character
add_compile_options($<$<Fortran_COMPILER_VERSION:9.3.0>:-Wno-maybe-uninitialized>)

check_fortran_compiler_flag(-fallow-argument-mismatch allow_mismatch_args)
if(allow_mismatch_args)
  set(gfortran_opts $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
endif()
