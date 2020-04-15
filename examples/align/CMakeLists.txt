add_executable(align_ex1 align_ex1_msa.cpp)
target_link_libraries(align_ex1 shared)

add_executable(align_ex2 align_ex2_libsimdpp.cpp)
target_link_libraries(align_ex2 shared)

add_executable(align_ex3 align_ex3_tests.cpp)
target_link_libraries(align_ex3 shared)

add_executable(align_ex4 align_ex4_input_sequences.cpp)
target_link_libraries(align_ex4 shared)

add_executable(align_ex5 align_ex5_graphs.cpp)
target_link_libraries(align_ex5 shared)

add_executable(align_ex6 align_ex6_variant_break_down.cpp)
target_link_libraries(align_ex6 shared)


# Search for clang-tidy
find_program(
    CLANG_TIDY_EXE
    NAMES "clang-tidy"
    DOC "Path to clang-tidy executable"
)

if(CLANG_TIDY_EXE)
  message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
  set(TIDY_CHECKS "clang-analyzer-*,clang-diagnostic-*,modernize-*,-clang-analyzer-alpha.*")
  set_target_properties(align_ex2 PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE};-checks=${TIDY_CHECKS}")
else()
  message(STATUS "clang-tidy not found. Code diagnostics and statis analyzes will be skipped.")
endif()
