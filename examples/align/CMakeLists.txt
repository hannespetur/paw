add_executable(align_ex1 ex1_msa.cpp)
target_link_libraries(align_ex1 shared)

add_executable(align_ex2 ex2_libsimdpp.cpp)
target_link_libraries(align_ex2 shared)

add_executable(align_ex3 ex3_benchmarks.cpp)
target_link_libraries(align_ex3 shared)

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
