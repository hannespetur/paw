add_executable(align_ex1 ex1_msa.cpp)
target_link_libraries(align_ex1 shared)

add_executable(align_ex2 ex2_libsimdpp.cpp)
target_link_libraries(align_ex2 shared)
