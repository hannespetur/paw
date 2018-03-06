//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================


#include <numeric>
#include <vector>

#include <boost/simd/algorithm/find.hpp>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_find()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(2 * N + 1);
  std::iota(values.begin(), values.end(), T(1));

  SECTION("t1")
  {
    auto f1 = std::find(values.begin(), values.end(), T(N));
    auto f2 = boost::simd::find(values.data(), values.data() + 2 * N + 1, T(N));
    REQUIRE(*f1 == *f2);
    auto f3 = std::find(values.begin(), values.end(), T(104));
    auto f4 = boost::simd::find(values.data(), values.data() + 2 * N + 1, T(104));
    REQUIRE(*f3 == *f4);
  }

  SECTION("t2")
  {
    auto f1 = std::find(values.begin(), values.end(), T(0));
    auto f2 = boost::simd::find(values.data(), values.data() + 2 * N + 1, T(0));
    REQUIRE(*f1 == *f2);
    auto f3 = std::find(values.begin(), values.end(), T(2*N));
    auto f4 = boost::simd::find(values.data(), values.data() + 2 * N + 1, T(2*N));
    REQUIRE(*f3 == *f4);
  }
}

TEST_CASE("test test_find")
{
  TEST_NUMERIC_TYPES(test_find);
}
