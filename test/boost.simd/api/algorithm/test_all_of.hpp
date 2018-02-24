#include <numeric>
#include <vector>

#include <boost/simd/algorithm/all_of.hpp>
#include <boost/align/aligned_allocator.hpp>

#include "../../../include/catch.hpp"

using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_all_of()
{
  static const int N = pack<T>::static_size;

  {
    std::vector<T, aligned_allocator<T, pack<T>::alignment>> values(static_cast<T>(2 * N));
    std::iota(values.begin(), values.end(), T(1));

    auto ab = values.data();
    auto ae = values.data() + values.size();

    REQUIRE(boost::simd::all_of(ab, ae));
    REQUIRE(boost::simd::all_of(ab + 1, ae));
    REQUIRE(boost::simd::all_of(ab + 1, ae));
    REQUIRE(boost::simd::all_of(ab, ae - 1));
    REQUIRE(boost::simd::all_of(ab, ae - 1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(static_cast<T>(2 * N));
    std::iota(values.begin(), values.end(),T(1));
    values[N] = T(0);  // 0 in aligned

    auto ab = values.data();
    auto ae = values.data()+values.size();

    REQUIRE_FALSE(boost::simd::all_of(ab, ae));
    REQUIRE_FALSE(boost::simd::all_of(ab+1, ae));
    REQUIRE_FALSE(boost::simd::all_of(ab, ae-1));
    REQUIRE_FALSE(boost::simd::all_of(ab+1, ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(static_cast<T>(2 * N));
    std::iota(values.begin(), values.end(),T(1));
    values[1] = T(0);  // 0 in prologue

    auto ab = values.data();
    auto ae = values.data()+values.size();

    REQUIRE_FALSE(boost::simd::all_of(ab,ae));
    REQUIRE_FALSE(boost::simd::all_of(ab+1,ae));
    REQUIRE_FALSE(boost::simd::all_of(ab,ae-1));
    REQUIRE_FALSE(boost::simd::all_of(ab+1,ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(static_cast<T>(2 * N));
    std::iota(values.begin(), values.end(),T(1));
    values[2*N-2] = T(0);// 0 in epilogue

    auto ab = values.data();
    auto ae = values.data()+values.size();

    REQUIRE_FALSE(boost::simd::all_of(ab,ae));
    REQUIRE_FALSE(boost::simd::all_of(ab+1,ae));
    REQUIRE_FALSE(boost::simd::all_of(ab,ae-1));
    REQUIRE_FALSE(boost::simd::all_of(ab+1,ae-1));
  }
}
