//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/algorithm/reduce.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <numeric>
#include <vector>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"


using namespace boost::simd;
using namespace boost::alignment;


struct fake_sum
{
  template<typename T> T operator()(T const& a, T const& e) { return a + e; }
};


template<typename T>
void
test_reduce()
{
  static const int N = pack<T>::static_size;

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(2*N);
    std::iota(values.begin(), values.end(),T(0));

    auto ab = values.data();
    auto ae = values.data()+values.size();

    // All aligned
    REQUIRE(std::accumulate(values.begin(), values.end(), T(3)) ==
            boost::simd::reduce(ab,ae, T(3))
            );

    // prologue + aligned
    REQUIRE(std::accumulate(values.begin() + 1, values.end(), T(3)) ==
            boost::simd::reduce(ab + 1, ae, T(3))
            );

    // aligned + epilogue
    REQUIRE(std::accumulate(values.begin(), values.end() - 1, T(3)) ==
            boost::simd::reduce(ab, ae - 1, T(3))
            );

    // prologue + epilogue
    REQUIRE(std::accumulate(values.begin() + 1, values.end() - 1, T(3)) ==
            boost::simd::reduce(ab, ae - 1, T(3))
            );
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(2 * N);
    std::iota(values.begin(), values.end(),T(0));

    REQUIRE(std::accumulate(values.begin(), values.end(), T(3), fake_sum{}) ==
            boost::simd::reduce(values.data(), values.data() + values.size(), T(3), fake_sum{})
            );
  }
}

TEST_CASE("test test_reduce")
{
  TEST_NUMERIC_TYPES(test_reduce);
}
