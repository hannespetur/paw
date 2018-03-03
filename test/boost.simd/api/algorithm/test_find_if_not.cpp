//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#pragma once

#include <numeric>
#include <vector>

#include <boost/simd/algorithm/find_if_not.hpp>

#include "../../../include/catch.hpp"
#include "help_structures.hpp"


using namespace boost::simd;


template<typename T>
void
test_find_if_not()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(2*N+1);
  std::iota(values.begin(), values.end(), T(1));

  {
    auto f1 = std::find_if_not(values.begin(), values.end(), [](T e) { return e < N; });
    auto f2 = boost::simd::find_if_not(values.data(), values.data()+2*N+1, f_lt<T>(N));
    REQUIRE(*f1 == *f2);
    auto f3 = std::find_if_not(values.begin(), values.end(), [](T e) { return e < 1024; });
    auto f4 = boost::simd::find_if_not(values.data(), values.data()+2*N+1, f_lt<T>(104));
    REQUIRE(*f3 == *f4);
  }

  {
    auto f1 = std::find_if_not(values.begin(), values.end(), [](T e) { return e < 0; });
    auto f2 = boost::simd::find_if_not(values.data(), values.data()+2*N+1, f_lt<T>(0));
    REQUIRE(*f1 == *f2);
    auto f3 = std::find_if_not(values.begin(), values.end(), [](T e) { return e < 2*N; });
    auto f4 = boost::simd::find_if_not(values.data(), values.data()+2*N+1, f_lt<T>(2*N));
    REQUIRE(*f3 == *f4);
  }
}
