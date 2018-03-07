#include <numeric>
#include <vector>

#include <boost/simd/algorithm/any_of.hpp>
#include <boost/align/aligned_allocator.hpp>

#include "../test_numeric_types.hpp"
#include <catch.hpp>


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_any_of()
{
  static const int N = pack<T>::static_size;

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(static_cast<T>(2 * N), T(0));

    auto ab = values.data();
    auto ae = values.data()+values.size();

    REQUIRE_FALSE(boost::simd::any_of(ab,ae));
    REQUIRE_FALSE(boost::simd::any_of(ab+1,ae));
    REQUIRE_FALSE(boost::simd::any_of(ab,ae-1));
    REQUIRE_FALSE(boost::simd::any_of(ab+1,ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(static_cast<T>(2 * N));
    std::iota(values.begin(), values.end(),T(1));
    values[N] = T(1); // 1 in aligned

    auto ab = values.data();
    auto ae = values.data()+values.size();

    REQUIRE(boost::simd::any_of(ab,ae));
    REQUIRE(boost::simd::any_of(ab+1,ae));
    REQUIRE(boost::simd::any_of(ab,ae-1));
    REQUIRE(boost::simd::any_of(ab+1,ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(static_cast<T>(2 * N));
    std::iota(values.begin(), values.end(),T(1));
    values[1] = T(1);// 1 in prologue

    auto ab = values.data();
    auto ae = values.data()+values.size();

    REQUIRE(boost::simd::any_of(ab,ae));
    REQUIRE(boost::simd::any_of(ab+1,ae));
    REQUIRE(boost::simd::any_of(ab,ae-1));
    REQUIRE(boost::simd::any_of(ab+1,ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(static_cast<T>(2 * N));
    std::iota(values.begin(), values.end(),T(1));
    values[2*N-2] = T(1);// 1 in epilogue

    auto ab = values.data();
    auto ae = values.data()+values.size();

    REQUIRE(boost::simd::any_of(ab,ae));
    REQUIRE(boost::simd::any_of(ab+1,ae));
    REQUIRE(boost::simd::any_of(ab,ae-1));
    REQUIRE(boost::simd::any_of(ab+1,ae-1));
  }
}

TEST_CASE("test any_of")
{
  TEST_NUMERIC_TYPES(test_any_of);
}

