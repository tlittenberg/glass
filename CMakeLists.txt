# glass

# This minimum version is mostly set in order to get a newer version
# of the FindMPI check.  Note that you can easily install a newer cmake version
# using conda or pip.
cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)

foreach(policy
    CMP0048
    CMP0074
    CMP0077
  )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

# Set compiler
set(CMAKE_C_COMPILER mpicc)

project(glass VERSION 2.0.0 LANGUAGES C)

# External packages
find_package(LAPACK REQUIRED)
find_package(OpenBLAS REQUIRED)
find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED)
find_package(MPI REQUIRED)
#find_package(GSL REQUIRED)

# set cflags
if(APPLE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xpreprocessor")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -Wall -std=gnu99")

# Set git version header file
message(STATUS "Resolving GIT Version")

find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE _build_version
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "GIT hash: ${_build_version}")
else()
  message(STATUS "GIT not found")
endif()


# The version number.
set (GITVERSION "\"${_build_version}\"")

configure_file (
  "${PROJECT_SOURCE_DIR}/gitversion.h.in"
  "${PROJECT_SOURCE_DIR}/utils/src/gitversion.h"
  )

install (FILES "${PROJECT_SOURCE_DIR}/utils/src/gitversion.h" DESTINATION include/glass)

message(STATUS "Resolving GIT Version")

# The version number.
set (GITVERSION "\"${_build_version}\"")

# Internal products
add_subdirectory(external)
add_subdirectory(utils)
add_subdirectory(ucb)
add_subdirectory(mbh)
add_subdirectory(noise)

# Merge into a single library
add_library(glass
        $<TARGET_OBJECTS:glass_external>
        $<TARGET_OBJECTS:glass_utils>
        $<TARGET_OBJECTS:glass_ucb>
        $<TARGET_OBJECTS:glass_mbh>
        $<TARGET_OBJECTS:glass_noise>
)
install(TARGETS glass DESTINATION lib)

# Applications
add_subdirectory(apps)
add_subdirectory(globalfit)
