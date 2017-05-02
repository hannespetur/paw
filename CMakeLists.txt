cmake_minimum_required (VERSION 3.2)

project (paw)

set (CMAKE_CXX_STANDARD 11)


if(NOT CMAKE_BUILD_TYPE STREQUAL Debug)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose: Debug or Release." FORCE)
endif()

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

#set (CMAKE_BUILD_TYPE CACHE RELEASE STRING "Build type")

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  message(STATUS "Building in debug mode")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG -pg -W -Wall -Wextra -Wfatal-errors -pedantic")
else ()
  message(STATUS "Building in release mode")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
endif ()

include_directories(include)

# Examples
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

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
endif()

# Testing
enable_testing()
add_subdirectory(test)
