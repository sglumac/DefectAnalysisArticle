# -----------------------------------------------------------------
# $Revision$
# $Date$
# -----------------------------------------------------------------
# Programmer: Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# CMakeLists.txt for  serial examples
#
# This file is generated from a template using  various variables
# set at configuration time. It can be used as a template for
# other user CMakeLists configuration files.
#
# Note: if the solver was successfully configured with Blas/Lapack
# support, the Blas/Lapack libraries are specified through the 
# variable LAPACK_LIBRARIES. Otherwise, this variable should contain
# an empty string or LAPACK_LIBRARIES-NOTFOUND. If this variable 
# contains a valid library entry, we add it to each example traget
# whether they use the Lapack module or not. This is done in order
# to address the case in which the SUNDIALS libraries are shared 
# objects. In that case, the solver library references Lapack 
# symbols which must be always resolved by linking against the
# Blas/Lapack libraries. If only static SUNDIALS libraries have 
# been built, it is not required to link the Blas/Lapack libraries
# for examples that do not use that module...
#
# -----------------------------------------------------------------

cmake_minimum_required(VERSION 2.8)

# Specify project name
PROJECT(_serial_examples C)

MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)

# Set the names of the examples to be built
SET(examples  test_nvector_serial )
LIST(REMOVE_DUPLICATES examples)

# Set names of examples specific dependency source files
SET(examples_dependencies  test_nvector sundials_nvector)

# Specify path to SUNDIALS header files
SET(SUNDIALS_INC_DIR
  C:/Users/slave/Desktop/Research Repositories/Co-simulation defects/sundials/include
  CACHE STRING
  "Location of SUNDIALS header files")

# Add path to SUNDIALS header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries 
SET(SUNDIALS_LIB_DIR C:/Users/slave/Desktop/Research Repositories/Co-simulation defects/sundials/lib)

# Find the SUNDIALS solver's library
FIND_LIBRARY(SUNDIALS_SOLVER_LIB
  sundials_nvecserial ${SUNDIALS_LIB_DIR}
  DOC " library")

# Find the NVECTOR library
FIND_LIBRARY(SUNDIALS_NVEC_LIB
  sundials_nvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR library")

# Set additional libraries
SET(SUNDIALS_EXTRA_LIB  CACHE STRING "Additional libraries")

# Set other libraries libraries
SET(LAPACK_LIBRARIES  CACHE STRING "Lapack libraries")
SET(SUPERLUMT_LIBRARIES  CACHE STRING "SuperLUMT libraries")
SET(KLU_LIBRARIES  CACHE STRING "KLU libraries")
SET(PETSC_LIBRARIES  CACHE STRING "PETSC libraries")

# List of all libraries
SET(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB} ${SUNDIALS_EXTRA_LIB})
IF(LAPACK_LIBRARIES)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_LIBRARIES)

IF(SUPERLUMT_LIBRARIES)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_LIBRARIES)

IF(KLU_LIBRARIES)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_LIBRARIES)

IF(PETSC_LIBRARIES)
  LIST(APPEND SUNDIALS_LIBS ${PETSC_LIBRARIES})
ENDIF(PETSC_LIBRARIES)

# Build each example one by one
FOREACH(example ${examples})
  ADD_EXECUTABLE(${example} ${example}.c ${examples_dependencies})
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
ENDFOREACH(example ${examples})
 
