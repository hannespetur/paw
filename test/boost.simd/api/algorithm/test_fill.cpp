//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================

#include <numeric>
#include <vector>

#include <boost/simd/algorithm/fill.hpp>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"

using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_fill()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values1(2*N+1), values2(2*N+1);;
  std::fill(values1.begin(), values1.end(),T(1));
  boost::simd::fill(values2.data(), values2.data()+values2.size(),T(1));

  REQUIRE(std::equal(values1.begin(), values1.end(), values2.begin()));
}

TEST_CASE("test test_fill")
{
  TEST_NUMERIC_TYPES(test_fill);
}
