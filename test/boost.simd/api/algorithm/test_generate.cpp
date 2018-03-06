//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/algorithm/generate.hpp>
#include <boost/simd/function/enumerate.hpp>
#include <boost/simd/pack.hpp>
#include <numeric>
#include <vector>

#include "../common.hpp"
#include "../../../include/catch.hpp"
#include "help_structures.hpp"


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_generate()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(2*N+1), ref(2*N+1);
  std::generate(ref.begin(), ref.end(), gstd());
  boost::simd::generate(values.data(), values.data()+values.size(), g());
  REQUIRE(values == ref);
  std::generate(ref.begin(), ref.end(), gstd(2, 3));
  boost::simd::generate(values.data(), values.data()+values.size(), g(2, 3));
  REQUIRE(values == ref);
}

TEST_CASE("test test_generate")
{
  TEST_NUMERIC_TYPES(test_generate);
}
