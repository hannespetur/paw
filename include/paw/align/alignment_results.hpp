#pragma once

#include <paw/align/libsimdpp_backtracker.hpp>

#include <simdpp/simd.h>


namespace paw
{

template <typename Tuint>
struct AlignmentResults
{
  long score{0};
  SIMDPP_ARCH_NAMESPACE::Backtrack<Tuint> mB;

  AlignmentResults() = default;
};

} // namespace paw