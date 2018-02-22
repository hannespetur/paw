//==================================================================================================
/**
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
**/
//==================================================================================================
#include <boost/simd/function/unary_minus.hpp>
#include <boost/simd/pack.hpp>
#include <simd_test.hpp>
#include <boost/simd/constant/inf.hpp>
#include <boost/simd/constant/minf.hpp>
#include <boost/simd/constant/mone.hpp>
#include <boost/simd/constant/nan.hpp>
#include <boost/simd/constant/one.hpp>
#include <boost/simd/constant/zero.hpp>
#include <boost/simd/constant/mzero.hpp>
#include <boost/simd/constant/two.hpp>
#
namespace bs = boost::simd;

template <typename T, std::size_t N, typename Env>
void test(Env& runtime)
{
  using p_t = bs::pack<T, N>;

  T a1[N], b[N];

  for(std::size_t i = 0; i < N; ++i)
  {
    a1[i] = (i%2) ? T(i) : T(2*i);
    b[i] = bs::unary_minus(a1[i]) ;
  }

  p_t aa1(&a1[0], &a1[0]+N);
  p_t bb (&b[0], &b[0]+N);

  STF_EQUAL(bs::unary_minus(aa1), bb);
  STF_EQUAL(-aa1, bb);
}

STF_CASE_TPL("Check unary_minus on pack" , STF_SIGNED_NUMERIC_TYPES)
{
  static const std::size_t N = bs::pack<T>::static_size;

  test<T, N>(runtime);
  test<T, N/2>(runtime);
  test<T, N*2>(runtime);
}


STF_CASE_TPL (" unary_minus real",  STF_IEEE_TYPES)
{
  namespace bs = boost::simd;
  namespace bd = boost::dispatch;
  using bs::unary_minus;
  using p_t = bs::pack<T>;

  using r_t = decltype(unary_minus(p_t()));

  // return type conformity test
  STF_TYPE_IS(r_t, p_t);

  // specific values tests
  STF_EQUAL(unary_minus(bs::Inf<p_t>()), bs::Minf<r_t>());
  STF_EQUAL(unary_minus(bs::Minf<p_t>()), bs::Inf<r_t>());
  STF_IEEE_EQUAL(unary_minus(bs::Nan<p_t>()), bs::Nan<r_t>());
  STF_EQUAL(unary_minus(bs::One<p_t>()), bs::Mone<r_t>());
  STF_EQUAL(unary_minus(bs::Zero<p_t>()), bs::Mzero<r_t>());
} // end of test for floating_

STF_CASE_TPL (" unary_minus signed_int",  STF_SIGNED_INTEGRAL_TYPES)
{
  namespace bs = boost::simd;
  namespace bd = boost::dispatch;
  using bs::unary_minus;
  using p_t = bs::pack<T>;
  using r_t = decltype(unary_minus(p_t()));

  // return type conformity test
  STF_TYPE_IS(r_t, p_t);

  // specific values tests
  STF_EQUAL(unary_minus(bs::Mone<p_t>()), bs::One<r_t>());
  STF_EQUAL(unary_minus(bs::One<p_t>()), bs::Mone<r_t>());
  STF_EQUAL(unary_minus(bs::Two<p_t>()), -bs::Two<r_t>());
  STF_EQUAL(unary_minus(bs::Zero<p_t>()), bs::Zero<r_t>());
} // end of test for signed_int_STF_CASE("unary_minus TO DO")

