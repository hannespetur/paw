cmake_minimum_required (VERSION 3.2)

project (pawparser)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_BUILD_TYPE CACHE Release STRING "Build type")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -fmax-errors=1")

include_directories(include)

enable_testing()
add_subdirectory(test)
