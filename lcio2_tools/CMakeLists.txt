cmake_minimum_required(VERSION 3.5)
project(LCIO2TOOLS CXX)

set(LCIO2TOOLS_MAJOR_VERSION 0)
set(LCIO2TOOLS_MINOR_VERSION 5)
set(LCIO2TOOLS_PATCH_VERSION 1)
set(LCIO2TOOLS_VERSION
  ${LCIO2TOOLS_MAJOR_VERSION}.${LCIO2TOOLS_MINOR_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

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

find_package(LCIO REQUIRED)

find_package(LCIO2 REQUIRED)
message(STATUS "Found LCIO2: ${LCIO2_DIR}")

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
  "${PROJECT_SOURCE_DIR}/include"
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_BINARY_DIR}/lcio"
  "${LCIO_INCLUDE_DIRS}"
  "${ROOT_INCLUDE_DIRS}"
  "${podio_INCLUDE_DIRS}"
  )

# --------------------------------------------------------------------------

file(GLOB sources src/lcio2utils.cc )
file(GLOB headers include/lcio2utils.h )

add_library(lcio2util SHARED ${sources} ${headers})
target_link_libraries(lcio2util lcio2 podio ${ROOT_LIBRARIES} ${LCIO_LIBRARIES})
set_target_properties(lcio2util PROPERTIES PUBLIC_HEADER "${headers}")

#REFLEX_GENERATE_DICTIONARY(lcio2 ${headers} SELECTION ${CMAKE_BINARY_DIR}/src/selection.xml )
#add_library(lcio2Dict SHARED lcio2.cxx)
#add_dependencies(lcio2Dict lcio2-dictgen)
#target_link_libraries(lcio2Dict lcio2 podio ${ROOT_LIBRARIES})
#
#set_target_properties(lcio2util PROPERTIES
#  PUBLIC_HEADER "${headers}")

install(TARGETS lcio2util
  # IMPORTANT: Add the datamodel library to the "export-set"
  #EXPORT fccedmTargets
  #RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/lcio2"
  COMPONENT dev)

#install(FILES
#  "${PROJECT_BINARY_DIR}/lcio2Dict.rootmap"
#  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT dev)

#if (${ROOT_VERSION} GREATER 6)
#  install(FILES
#      "${PROJECT_BINARY_DIR}/lcio2_rdict.pcm"
#      DESTINATION "${INSTALL_LIB_DIR}" COMPONENT dev)
#endif()

# --------------------------------------------------------------------------

add_executable(legacy2lcio2 src/legacy2lcio2.cc)
target_link_libraries(legacy2lcio2 lcio2util lcio2 podio ${ROOT_LIBRARIES} ${LCIO_LIBRARIES})
add_dependencies(legacy2lcio2 lcio2util)

add_executable(lcio2legacy src/lcio2legacy.cc)
target_link_libraries(lcio2legacy lcio2util lcio2 podio ${ROOT_LIBRARIES} ${LCIO_LIBRARIES})
add_dependencies(lcio2legacy lcio2util)


add_executable(lcio2-edm-read src/read.cc)
target_link_libraries(lcio2-edm-read lcio2util lcio2 podio ${ROOT_LIBRARIES} ${LCIO_LIBRARIES})
add_dependencies(lcio2-edm-read lcio2util)

install(TARGETS lcio2-edm-read legacy2lcio2 lcio2legacy
  #EXPORT fccedmTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)

# --------------------------------------------------------------------------

