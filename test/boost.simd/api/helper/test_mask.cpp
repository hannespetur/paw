//==================================================================================================
/*
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/mask.hpp>
#include <boost/simd/detail/dispatch/hierarchy_of.hpp>

#include <catch.hpp>


namespace bs = boost::simd;
namespace bd = boost::dispatch;


TEST_CASE( "Check masked pointer interface" )
{
  double x {};
  std::uint64_t const const_x{};

  auto ax = bs::mask(&x, true, 69.);
  auto acx = bs::mask(&const_x, false, 42ul);

  REQUIRE(ax.get() == &x);
  REQUIRE(acx.get() == &const_x);
  REQUIRE(ax.mask());
  REQUIRE_FALSE(acx.mask());
  REQUIRE(ax.value() == 69.);
  REQUIRE(acx.value() == 42u);
}


TEST_CASE("hierarchy_of of masked pointer")
{
  float f{};
  char const c{};

  auto m1 = bs::mask(&f, true, f+1);
  auto m2 = bs::mask(&f, true);
  auto m3 = bs::mask(&c, false, c+1);
  auto m4 = bs::mask(&c, false);

  using masked_ptr_t        = decltype(m1);
  using zero_masked_ptr_t   = decltype(m2);
  using c_masked_ptr_t      = decltype(m3);
  using zero_c_masked_ptr_t = decltype(m4);

  {
    using T = bd::hierarchy_of_t<masked_ptr_t>;
    using U = bd::hierarchy_of_t<float, masked_ptr_t>;
    bool constexpr b = std::is_same<T, bd::masked_pointer_<U, std::false_type> >::value;
    REQUIRE(b);
  }

  {
    using T = bd::hierarchy_of_t<zero_masked_ptr_t>;
    using U = bd::hierarchy_of_t<float, zero_masked_ptr_t>;
    bool constexpr b = std::is_same<T, bd::masked_pointer_<U, std::true_type> >::value;
    REQUIRE(b);
  }

  {
    using T = bd::hierarchy_of_t<c_masked_ptr_t>;
    using U = bd::hierarchy_of_t<char const, c_masked_ptr_t>;
    bool constexpr b = std::is_same<T, bd::masked_pointer_<U, std::false_type> >::value;
    REQUIRE(b);
  }

  {
    using T = bd::hierarchy_of_t<zero_c_masked_ptr_t>;
    using U = bd::hierarchy_of_t<char const, zero_c_masked_ptr_t>;
    bool constexpr b = std::is_same<T, bd::masked_pointer_<U, std::true_type> >::value;
    REQUIRE(b);
  }
}
