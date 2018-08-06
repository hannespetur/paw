#pragma once

#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>

#include <simdpp/simd.h>


namespace paw
{

template <typename Tuint>
struct AlignmentResults
{
  long score{0};
  SIMDPP_ARCH_NAMESPACE::Backtrack<Tuint> mB;
  typename T<Tuint>::vec_pack last_vH;
  typename T<Tuint>::vec_pack last_vF;

  AlignmentResults() = default;
};

} // namespace paw