cmake_minimum_required(VERSION 2.8)

project(SPIN_BOSON Fortran)

set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_Fortran_FLAGS   "-O0 ")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

set(TARGET_NAME   "spin-boson")
set(MOD_LIB        module)
set(IO_LIB         io)
set(COMMON_LIB     common)
set(GS_LIB         gs)
set(RT_LIB         rt)

file(GLOB SOURCES src/*.f90)

add_executable(${TARGET_NAME} ${SOURCES})
