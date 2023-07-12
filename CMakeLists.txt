cmake_minimum_required (VERSION 2.6)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
project (MDLovoFit Fortran)
add_executable(mdlovofit 
    src/common.f90
    src/lib.f90
    src/mdlovofit.f90
)

# Installation directive
install(TARGETS mdlovofit DESTINATION bin)



