//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/algorithm/max_element.hpp>
#include <boost/simd/function/is_greater.hpp>
#include <numeric>
#include <vector>

#include "../../../include/catch.hpp"

using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_max_element()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(2*N+1);
  {
    std::iota(values.begin(), values.end(), T(1));
    values[N] = T(1000);
    auto f1 = std::max_element(values.begin(), values.end());
    auto f2 = boost::simd::max_element(values.data(), values.data()+2*N+1);
    REQUIRE(*f1 == *f2);
  }

  {
    std::iota(values.begin(), values.end(), T(1));
    values[0] = T(1000);
    auto f1 = std::max_element(values.begin(), values.end());
    auto f2 = boost::simd::max_element(values.data(), values.data()+2*N+1);
    REQUIRE(*f1 == *f2);
  }

  {
    std::iota(values.begin(), values.end(), T(1));
    values[2*N] = T(0);
    auto f1 = std::max_element(values.begin(), values.end());
    auto f2 = boost::simd::max_element(values.data(), values.data()+2*N+1);
    REQUIRE(*f1 == *f2);
  }

  {
    std::iota(values.begin(), values.end(), T(1));
    values[N] = T(0);
    auto f1 = std::max_element(values.begin(), values.end(), bs::is_greater);
    auto f2 = boost::simd::max_element(values.data(), values.data()+2*N+1, bs::is_greater);
    REQUIRE(*f1 == *f2);
  }

  {
    std::iota(values.begin(), values.end(), T(1));
    values[0] = T(0);
    auto f1 = std::max_element(values.begin(), values.end(), bs::is_greater);
    auto f2 = boost::simd::max_element(values.data(), values.data()+2*N+1, bs::is_greater);
    REQUIRE(*f1 == *f2);
  }

  {
    std::iota(values.begin(), values.end(), T(1));
    values[2*N] = T(0);
    auto f1 = std::max_element(values.begin(), values.end(), bs::is_greater);
    auto f2 = boost::simd::max_element(values.data(), values.data()+2*N+1, bs::is_greater);
    REQUIRE(*f1 == *f2);
  }
}
