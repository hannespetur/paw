cmake_minimum_required (VERSION 2.8.12)

project (paw)

## Set to release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE is set to ${CMAKE_BUILD_TYPE}")

## Optional dependencies
# Boost
find_package(Boost COMPONENTS filesystem system)

## Setup compile flags
set (CMAKE_CXX_FLAGS_COMMON "-g -std=gnu++11 -Wall -Wextra -Wfatal-errors")

if(NOT PORTABLE)
  set (CMAKE_CXX_FLAGS_COMMON "${CMAKE_CXX_FLAGS_COMMON} -march=native")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_COMMON}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_COMMON}")

## Build libraries in /lib and binaries in /bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Include paw
include_directories("${PROJECT_SOURCE_DIR}/include")

## Static and shared libraries
# Note: You can make a static library with 'make static' and a shared library with 'make shared'. These libraries will have ready all implementations (of non-templated functions). If you use you should not define IMPLEMENT in your any cpp file as you would otherwise need)
add_library(static STATIC EXCLUDE_FROM_ALL paw.cpp)
SET_TARGET_PROPERTIES(static PROPERTIES OUTPUT_NAME paw)

add_library(shared SHARED paw.cpp)
SET_TARGET_PROPERTIES(shared PROPERTIES OUTPUT_NAME paw)
set_property(TARGET shared PROPERTY POSITION_INDEPENDENT_CODE ON)

## Examples
add_subdirectory(examples)

# Documentation
## Adds a target to generate API documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()

  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

  #install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
endif()

## Testing
enable_testing()
add_subdirectory(test)
