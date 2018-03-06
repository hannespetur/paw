//==================================================================================================
/*
  Copyright 2017 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/algorithm/none_of.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <numeric>
#include <vector>

#include "../common.hpp"
#include "../../../include/catch.hpp"


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_none_of()
{
  static const int N = pack<T>::static_size;

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(2*N, T(0));

    auto ab = values.data();
    auto ae = values.data() + values.size();

    // All aligned
    REQUIRE(boost::simd::none_of(ab,ae));

    // prologue + aligned
    REQUIRE(boost::simd::none_of(ab+1,ae));

    // aligned + epilogue
    REQUIRE(boost::simd::none_of(ab,ae-1));

    // prologue + epilogue
    REQUIRE(boost::simd::none_of(ab+1,ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(2*N);
    std::iota(values.begin(), values.end(),T(1));
    values[N] = T(1); // 1 in aligned

    auto ab = values.data();
    auto ae = values.data()+values.size();

    // All aligned
    REQUIRE_FALSE(boost::simd::none_of(ab,ae));

    // prologue + aligned
    REQUIRE_FALSE(boost::simd::none_of(ab+1,ae));

    // aligned + epilogue
    REQUIRE_FALSE(boost::simd::none_of(ab,ae-1));

    // prologue + epilogue
    REQUIRE_FALSE(boost::simd::none_of(ab+1,ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(2*N);
    std::iota(values.begin(), values.end(),T(1));
    values[1] = T(1);// 1 in prologue

    auto ab = values.data();
    auto ae = values.data()+values.size();

    // All aligned
    REQUIRE_FALSE (boost::simd::none_of(ab,ae));

    // prologue + aligned
    REQUIRE_FALSE (boost::simd::none_of(ab+1,ae));

    // aligned + epilogue
    REQUIRE_FALSE(boost::simd::none_of(ab,ae-1));

    // prologue + epilogue
    REQUIRE_FALSE(boost::simd::none_of(ab+1,ae-1));
  }

  {
    std::vector<T,aligned_allocator<T,pack<T>::alignment>> values(2*N);
    std::iota(values.begin(), values.end(),T(1));
    values[2*N-2] = T(1);// 1 in epilogue

    auto ab = values.data();
    auto ae = values.data()+values.size();

    // All aligned
    REQUIRE_FALSE(boost::simd::none_of(ab,ae));

    // prologue + aligned
    REQUIRE_FALSE(boost::simd::none_of(ab+1,ae));

    // aligned + epilogue
    REQUIRE_FALSE(boost::simd::none_of(ab,ae-1));

    // prologue + epilogue
    REQUIRE_FALSE(boost::simd::none_of(ab+1,ae-1));
  }
}

TEST_CASE("test test_none_of")
{
  TEST_NUMERIC_TYPES(test_none_of);
}
