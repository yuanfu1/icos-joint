cmake_minimum_required(VERSION 3.1)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(icos)

# Include the common template:
IF (IS_ABSOLUTE $ENV{CTEMPLATES})
  MESSAGE(STATUS "CTEMPLATES path is set!")
ELSE ($ENV{CTEMPLATES} LESS 1)
  MESSAGE(FATAL_ERROR "Env CTEMPLATES is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{CTEMPLATES})
include($ENV{CTEMPLATES}/cmakeCommon.txt)

#-------------------------------------------------------------------------------
# Source codes to this build:
include (FilesSelected.txt)

#-------------------------------------------------------------------------------

ADD_LIBRARY(linpack STATIC IMPORTED)
ADD_LIBRARY(blas STATIC IMPORTED)
SET_TARGET_PROPERTIES(linpack PROPERTIES IMPORTED_LOCATION $ENV{HOME}/developments/sftware/netlib/LINPACK/lib/liblinpack.a)
SET_TARGET_PROPERTIES(blas PROPERTIES IMPORTED_LOCATION $ENV{HOME}/developments/sftware/netlib/LINPACK/lib/liblas.a)

# Compiler options:
IF (Fortran_COMPILER_NAME MATCHES "pgf")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-D_REAL4_")
  SET (CMAKE_Fortran_FLAGS_DEBUG "-D_REAL4_")
ENDIF (Fortran_COMPILER_NAME MATCHES "pgf")
IF (Fortran_COMPILER_NAME MATCHES "gfortran")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-fbounds-check -O3 -freal-4-real-8 -D USE_EDGE_STCL -D DIAG")
  SET (CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -DDEBUG -freal-4-real-8 -D USE_EDGE_STCL -D DIAG")
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG}")
MESSAGE(STATUS "XXIE" ${CMAKE_Fortran_FLAGS_DEBUG})
ENDIF (Fortran_COMPILER_NAME MATCHES "gfortran")

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(icos STATIC ${icos_SRC})

# 2. Build multigrid generator executable:
SET (EXE1 "ml_itmesh.exe")

# 3. Add code to the executables:
ADD_EXECUTABLE(${EXE1} ml_itmesh.F90)

#-------------------------------------------------------------------------------
#  To link another lib like -l in cmake:
#-------------------------------------------------------------------------------
#  TARGET_LINK_LIBRARIES(target <list of libnames>)
TARGET_LINK_LIBRARIES(${EXE1} icos linpack blas)
