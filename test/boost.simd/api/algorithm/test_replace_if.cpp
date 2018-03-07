//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/algorithm/replace_if.hpp>
#include <boost/simd/pack.hpp>
#include <numeric>
#include <vector>

#include "../test_numeric_types.hpp"
#include <catch.hpp>


namespace bs = boost::simd;


template < typename T>
struct p
{
  using p_t = bs::pack<T>;

  p(const T & val)
    : val_(val)
    , pval_(val)
  {}

  T val_;
  p_t pval_;

  bool                  operator ()(T a)          const { return a < val_; }
  bs::as_logical_t<p_t> operator ()(const p_t& a) const { return a < pval_; }
};


template<typename T>
void
test_replace_if()
{
  static const int N = bs::pack<T>::static_size;

  std::vector<T> values(2 * N + 3), ref(2 * N + 3);
  std::iota(values.begin(), values.end(), T(1));
  std::iota(ref.begin(), ref.end(), T(1));

  SECTION("t1")
  {
    std::replace_if(ref.begin(), ref.end(), [](T e){ return e < N; }, T(0));
    boost::simd::replace_if(values.data(), values.data() + 2 * N + 3, p<T>(N), T(0));
    REQUIRE(values == ref);
  }

  SECTION("t2")
  {
    std::replace_if(ref.begin(), ref.end(), [](T e){ return e < 1; }, T(0));
    boost::simd::replace_if(values.data(), values.data() + 2 * N + 3, p<T>(1), T(0));
    REQUIRE(values == ref);
  }

  SECTION("t3")
  {
    std::replace_if(ref.begin(), ref.end(),  [](T e){ return e < 2*N-1; }, T(0));
    boost::simd::replace_if(values.data(), values.data() + 2 * N + 3, p<T>(2 * N - 1), T(0));
    REQUIRE(values == ref);
  }
}


TEST_CASE("test test_replace_if")
{
  TEST_NUMERIC_TYPES(test_replace_if);
}
