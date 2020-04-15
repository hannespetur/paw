add_library(catch SHARED main.cpp)
set_property(TARGET catch PROPERTY POSITION_INDEPENDENT_CODE ON)
list(APPEND PAW_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/test/include")

#add_subdirectory(align)
#add_subdirectory(boost.simd)
add_subdirectory(parser)
add_subdirectory(station)

# Use 'make check' to get verbose test results
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose DEPENDS test_pawparser test_pawstation test_align)
