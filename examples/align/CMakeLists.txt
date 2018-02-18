cmake_minimum_required (VERSION 2.8.11)

add_executable(align_ex1 ex1_msa.cpp)
target_link_libraries("align_ex1" shared ${PAW_LINK_LIBRARIES})
