cmake_minimum_required (VERSION 3.2)

project (paw CXX)

include(CheckCXXCompilerFlag)

# Include custom CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(SimdppMultiarch)


# Find which archs are compilable and which are runnable
simdpp_get_compilable_archs(COMPILABLE_ARCHS)
simdpp_get_runnable_archs(RUNNABLE_ARCHS)

message(STATUS "Compilable archs are: ${COMPILABLE_ARCHS}")
message(STATUS "Runnable archs are: ${RUNNABLE_ARCHS}")


set(PAW_VERSION_MAJOR 0)
set(PAW_VERSION_MINOR 1)
set(PAW_VERSION_PATCH 0)

# Get the current working git branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PAW_GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated SHA commit hash of the working git branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PAW_GIT_COMMIT_SHORT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest SHA commit hash of the working git branch
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PAW_GIT_COMMIT_LONG_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest SHA commit hash of the working git branch
execute_process(
  COMMAND git diff-files COMMAND wc -l
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PAW_GIT_NUM_DIRTY_FILES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

###
# OPTIONAL DEPENDENCIES
###
# Find Boost
find_package(Boost COMPONENTS filesystem system iostreams)

if (Boost_FOUND)
  message(STATUS "Including boost directory: ${Boost_INCLUDE_DIRS}")
  list(APPEND PAW_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})

  message(STATUS "Linking boost libraries: ${Boost_LIBRARIES}")
  list(APPEND PAW_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()


# Pass information of paw to the source
configure_file(${PROJECT_SOURCE_DIR}/include/paw/internal/_config.hpp
               ${PROJECT_BINARY_DIR}/include/paw/internal/config.hpp
               @ONLY)

# Include configured files
list(APPEND PAW_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/include)

## Set to release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE is set to ${CMAKE_BUILD_TYPE}")

## Find zlib
find_package(ZLIB REQUIRED)
list(APPEND PAW_LINK_LIBRARIES ${ZLIB_LIBRARIES})

# Find thread library (e.g. pthread on unix)
find_package (Threads)
list(APPEND PAW_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# Check which compiler flags are supported
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)

# Use highest CXX compiler flag available
if (COMPILER_SUPPORTS_CXX17)# AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Compilation fails on GCC 7.3 with c++17.
  set(CMAKE_CXX_FLAGS_COMMON "-std=c++17")
elseif (COMPILER_SUPPORTS_CXX14)
  set(CMAKE_CXX_FLAGS_COMMON "-std=c++14")
elseif (COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS_COMMON "-std=c++11")
else()
  message(FATAL_ERROR "Your compiler version does not seem to support C++11. Aborting.")
endif()

## Setup other compiler flags
set(CMAKE_CXX_FLAGS_COMMON "${CMAKE_CXX_FLAGS_COMMON} -Wall -Wextra -Wfatal-errors -m64 -march=x86-64 -mtune=generic")

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DPAW_BUILD ${CMAKE_CXX_FLAGS_COMMON} -pedantic -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -DPAW_BUILD ${CMAKE_CXX_FLAGS_COMMON} -DNDEBUG")

## Build libraries in /lib and binaries in /bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#simdpp_multiarch(GEN_ARCH_FILES src/align.cpp "NONE_NULL;X86_SSE2;X86_SSE3;X86_SSE4_1;X86_SSE4_1,X86_POPCNT_INSN;X86_AVX,X86_POPCNT_INSN;X86_AVX2,X86_POPCNT_INSN;X86_AVX512F,X86_POPCNT_INSN;X86_AVX512F,X86_POPCNT_INSN,X86_AVX512BW,X86_AVX512DQ,X86_AVX512VL") # ${RUNNABLE_ARCHS}, ${COMPILABLE_ARCHS}

if (FORCE_AVX)
  simdpp_multiarch(GEN_ARCH_FILES src/align.cpp "NONE_NULL;X86_SSE2;X86_SSE4_1;X86_SSE4_1,X86_POPCNT_INSN;X86_AVX,X86_POPCNT_INSN")
elseif (FORCE_AVX2)
  simdpp_multiarch(GEN_ARCH_FILES src/align.cpp "NONE_NULL;X86_SSE2;X86_SSE4_1;X86_SSE4_1,X86_POPCNT_INSN;X86_AVX,X86_POPCNT_INSN;X86_AVX2,X86_POPCNT_INSN")
elseif (FORCE_AVX512)
  simdpp_multiarch(GEN_ARCH_FILES src/align.cpp "NONE_NULL;X86_SSE2;X86_SSE4_1;X86_SSE4_1,X86_POPCNT_INSN;X86_AVX,X86_POPCNT_INSN;X86_AVX2,X86_POPCNT_INSN;X86_AVX512F,X86_POPCNT_INSN,X86_AVX512BW,X86_AVX512DQ,X86_AVX512VL")
else ()
  simdpp_multiarch(GEN_ARCH_FILES src/align.cpp ${COMPILABLE_ARCHS})
endif()

#add_executable(${PROJECT_NAME} src/main.cpp ${GEN_ARCH_FILES})

## Static and shared libraries
# Note: You can make a static library with 'make static' and a shared library with 'make shared'. These libraries will have ready all implementations (of non-templated functions). If you use you should not define IMPLEMENT in your any cpp file as you would otherwise need)
add_library(static STATIC EXCLUDE_FROM_ALL src/paw.cpp ${GEN_ARCH_FILES})
set_target_properties(static PROPERTIES
    OUTPUT_NAME paw
    )
target_link_libraries(static INTERFACE ${PAW_LINK_LIBRARIES})
target_include_directories(static SYSTEM PUBLIC ${PAW_INCLUDE_DIRS})
target_include_directories(static PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_library(shared SHARED src/paw.cpp ${GEN_ARCH_FILES})
set_target_properties(shared PROPERTIES OUTPUT_NAME paw POSITION_INDEPENDENT_CODE ON)
target_link_libraries(shared INTERFACE ${PAW_LINK_LIBRARIES})
target_include_directories(shared SYSTEM PUBLIC ${PAW_INCLUDE_DIRS})
target_include_directories(shared PUBLIC ${PROJECT_SOURCE_DIR}/include)


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

## Install
install(TARGETS shared LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include PATTERN include/paw/internal EXCLUDE)
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/paw/internal/ DESTINATION include/paw/internal)

## Testing
enable_testing()
add_subdirectory(test)
