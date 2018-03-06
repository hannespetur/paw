//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/algorithm/replace.hpp>
#include <numeric>
#include <vector>

#include "../common.hpp"
#include "../../../include/catch.hpp"


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_replace()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(2*N+3), ref(2*N+3);
  std::iota(values.begin(), values.end(), T(1));
  std::iota(ref.begin(), ref.end(), T(1));

  {
    std::replace(ref.begin(), ref.end(), T(1), T(0));
    boost::simd::replace(values.data(), values.data() + 2 * N + 3, T(1), T(0));
    REQUIRE(values == ref);
  }

  {
    std::replace(ref.begin(), ref.end(), T(2*N+2), T(0));
    boost::simd::replace(values.data(), values.data() + 2 * N + 3, T(2*N+2), T(0));
    REQUIRE(values == ref);
  }

  {
    std::replace(ref.begin(), ref.end(), T(2*N+2), T(0));
    boost::simd::replace(values.data(), values.data() + 2 * N + 3, T(2*N+2), T(0));
    REQUIRE(values == ref);
  }
}

TEST_CASE("test test_replace")
{
  TEST_NUMERIC_TYPES(test_replace);
}
