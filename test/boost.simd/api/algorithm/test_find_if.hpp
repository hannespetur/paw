//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================

#include <numeric>
#include <vector>

#include <boost/simd/algorithm/find_if.hpp>

#include "../../../include/catch.hpp"
#include "f_struct.hpp"


using namespace boost::simd;


template<typename T>
void
test_find_if()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(2*N+1);
  std::iota(values.begin(), values.end(), T(1));

  {
    auto f1 = std::find_if(values.begin(), values.end(), [](T e) { return e >= N; } );
    auto f2 = boost::simd::find_if(values.data(), values.data()+2*N+1, f_ge<T>(T(N)));
    REQUIRE(*f1 == *f2);
    auto f3 = std::find_if(values.begin(), values.end(), [](T e) { return e >= 104; } );
    auto f4 = boost::simd::find_if(values.data(), values.data()+2*N+1, f_ge<T>(T(104)) );
    REQUIRE(*f3 == *f4);
  }

  {
    auto f1 = std::find_if(values.begin(), values.end(), [](T e) { return e >= 0; } );
    auto f2 = boost::simd::find_if(values.data(), values.data()+2*N+1, f_ge<T>(T(0)) );
    REQUIRE(*f1 == *f2);
    auto f3 = std::find_if(values.begin(), values.end(), [](T e) { return e >= 2*N; });
    auto f4 = boost::simd::find_if(values.data(), values.data()+2*N+1, f_ge<T>(T(2*N)));
    REQUIRE(*f3 == *f4);
  }
}
