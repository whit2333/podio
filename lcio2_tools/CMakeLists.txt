cmake_minimum_required(VERSION 3.5)
project(LCIO2TOOLS CXX)

set(LCIO2TOOLS_MAJOR_VERSION 0)
set(LCIO2TOOLS_MINOR_VERSION 5)
set(LCIO2TOOLS_PATCH_VERSION 1)
set(LCIO2TOOLS_VERSION
  ${LCIO2TOOLS_MAJOR_VERSION}.${LCIO2TOOLS_MINOR_VERSION})

#execute_process(COMMAND podio ${CMAKE_SOURCE_DIR}/lcio-edm.yaml ${PROJECT_BINARY_DIR} lcio )

#--- Declare options -----------------------------------------------------------
#option(CREATE_DOC "Whether or not to create doxygen doc target.")

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")


if(APPLE)
  set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

set(DEF_INSTALL_CMAKE_DIR cmake)
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Set up C++14
set(CPP_STANDARD_FLAGS "-std=c++14")
if (${APPLE})
    set(CPP_STANDARD_FLAGS "-std=c++14\ -stdlib=libc++")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DDROP_CGAL ${CPP_STANDARD_FLAGS} -Wall -Wextra -Wpedantic -Wno-unused-variable -Wno-unused-parameter")


# Make sure we find the Find*.cmake functions distributed with this package
#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(LCIO2 REQUIRED)
message(STATUS "Found LCIO2: ${LCIO2_DIR}")
# Make sure the library is found.
# Not the case if LD_LIBRARY_PATH is wrong
#find_library(LCIO2LIB lcio2)

link_directories(${podio_LIBRARY_DIR})
# Make sure we find the Find*.cmake functions distributed with this package
#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(podio REQUIRED)
message(STATUS "Found podio: ${podio_DIR}")
# Make sure the library is found.
# Not the case if LD_LIBRARY_PATH is wrong
find_library(PODIOLIB podio PATHS $ENV{PODIO}/lib)
#if (NOT PODIOLIB)
#    message(FATAL_ERROR "libpodio.so(dylib) cannot be found dynamically. Make sure you have sourced PODIO init*.sh file to set up your environment to use PODIO")
#endif()
link_directories(${podio_LIBRARY_DIR})

#--- Declare ROOT dependency ---------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(ROOT REQUIRED COMPONENTS RIO Tree Physics)
include_directories(${ROOT_INCLUDE_DIR})
include(${ROOT_USE_FILE})

#find_package(dag REQUIRED HINTS $ENV{FCCDAG})
#if(${DAG_FOUND})
#  message(STATUS "Found dag: ${DAG_DIR}")
#  include_directories(
#          ${DAG_INCLUDE_DIRS}
#  )
#endif(${DAG_FOUND})

link_directories(${ROOT_LIBRARY_DIR})
#add_definitions(-Wpadded)
#--temporary fix of inconsistency in ROOT CMake macros
set(ROOT_genreflex_cmd ${ROOT_genreflex_CMD})

add_definitions(-Wno-unused-variable -Wno-unused-parameter)

message(
  "${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
  ${ROOT_INCLUDE_DIR}
  ${podio_INCLUDE_DIRS}"
  )

# set up include-directories
include_directories(
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_BINARY_DIR}/lcio"
  "${ROOT_INCLUDE_DIRS}"
  "${podio_INCLUDE_DIRS}"
  )


# adding testing capabilities
#include(CTest)

#--- target for Doxygen documentation ------------------------------------------
#jif(LCIO2TOOLS_documentation)
#j  include(cmake/LCIO2TOOLSDoxygen.cmake)
#jendif()

#configure_file(${CMAKE_SOURCE_DIR}/cmake/Datamodel.cmake
#  "${CMAKE_SOURCE_DIR}/datamodel/CMakeLists.txt" @ONLY)
# Add sub-directories
#add_subdirectory(datamodel)
#add_subdirectory(utilities)
#add_subdirectory(examples)
#add_subdirectory(tests)

#file(GLOB sources ${CMAKE_BINARY_DIR}/src/*.cc)
#file(GLOB headers ${CMAKE_BINARY_DIR}/lcio/*.h podio/PythonEventStore.h)
#
#add_library(LCIO2TOOLS SHARED ${sources} ${headers})
#target_link_libraries(LCIO2TOOLS podio)
#
#REFLEX_GENERATE_DICTIONARY(LCIO2TOOLS ${headers} SELECTION ${CMAKE_BINARY_DIR}/src/selection.xml )
#add_library(LCIO2TOOLSDict SHARED LCIO2TOOLS.cxx)
#add_dependencies(LCIO2TOOLSDict datamodel-dictgen)
#target_link_libraries(LCIO2TOOLSDict LCIO2TOOLS podio ${ROOT_LIBRARIES})
#
#set_target_properties(LCIO2TOOLS PROPERTIES
#  PUBLIC_HEADER "${headers}")
#
#install(TARGETS LCIO2TOOLS LCIO2TOOLSDict
#  # IMPORTANT: Add the datamodel library to the "export-set"
#  #EXPORT fccedmTargets
#  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
#  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
#  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/LCIO2TOOLS"
#  COMPONENT dev)
#
#install(FILES
#  "${PROJECT_BINARY_DIR}/LCIO2TOOLSDict.rootmap"
#  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT dev)
#
#if (${ROOT_VERSION} GREATER 6)
#  install(FILES
#      "${PROJECT_BINARY_DIR}/LCIO2TOOLS_rdict.pcm"
#      DESTINATION "${INSTALL_LIB_DIR}" COMPONENT dev)
#endif()
#
#
### The interesting stuff goes here
### ===============================
##
### Add all targets to the build-tree export set
##export(TARGETS datamodel utilities LCIO2TOOLS-write LCIO2TOOLS-read LCIO2TOOLS-simplewrite
##  FILE "${PROJECT_BINARY_DIR}/LCIO2TOOLSTargets.cmake")
##
### Export the package for use from the build-tree
### (this registers the build-tree with a global CMake-registry)
##export(PACKAGE LCIO2TOOLS)
##
### Create the LCIO2TOOLSConfig.cmake and LCIO2TOOLSConfigVersion files
##file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
##   "${INSTALL_INCLUDE_DIR}")
### ... for the build tree
##set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
#
#configure_file(cmake/LCIO2TOOLSConfig.cmake.in
#  "${PROJECT_BINARY_DIR}/LCIO2TOOLSConfig.cmake" @ONLY)
#
#configure_file(cmake/LCIO2TOOLSConfigVersion.cmake.in
#  "${PROJECT_BINARY_DIR}/LCIO2TOOLSConfigVersion.cmake" @ONLY)
#
## Install the LCIO2TOOLSConfig.cmake and LCIO2TOOLSConfigVersion.cmake
#install(FILES
#  "${PROJECT_BINARY_DIR}/LCIO2TOOLSConfig.cmake"
#  "${PROJECT_BINARY_DIR}/LCIO2TOOLSConfigVersion.cmake"
#  DESTINATION lib/LCIO2TOOLS)
##
### Install the export set for use with the install-tree
##install(EXPORT LCIO2TOOLSTargets DESTINATION
##  "${INSTALL_CMAKE_DIR}" COMPONENT dev)
