//==================================================================================================
/*
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/logical.hpp>

#include "../test_numeric_types.hpp"
#include <catch.hpp>

using namespace boost::simd;
namespace bs = boost::simd;

template<typename T>
void
test_logical()
{
  bs::logical<T> bool_true(true);
  bs::logical<T> bool_false(false);
  bs::logical<T> val_true(T(42));
  bs::logical<T> val_false(T(0));

  SECTION("Check logical<T> conversion from/to bool")
  {
    REQUIRE(bool(bool_true));
    REQUIRE(bool_true.value());
    REQUIRE_FALSE(bool(bool_false));
    REQUIRE_FALSE(bool_false.value());
    REQUIRE(bool(val_true));
    REQUIRE(val_true.value());
    REQUIRE_FALSE(bool(val_false));
    REQUIRE_FALSE(val_false.value());
  }

  SECTION("Check if and ?:")
  {
    if(val_false)
      REQUIRE(false); // logical<T> dont behave well in if()

    if((val_true ? false  : true))
      REQUIRE(false); // logical<T> dont behave well in ?:
  }

  SECTION("Check logical<T> comparison operators")
  {
    REQUIRE_FALSE((!bool_true).value());
    REQUIRE_FALSE((~bool_true).value());
    REQUIRE((!bool_false).value());
    REQUIRE((~bool_false).value());
    REQUIRE_FALSE((!val_true).value());
    REQUIRE_FALSE((~val_true).value());
    REQUIRE((!val_false).value());
    REQUIRE((~val_false).value());
  }

  SECTION("Check logical comparisons")
  {
    REQUIRE(bool_true  == val_true );
    REQUIRE(bool_false == val_false);
    REQUIRE(bool_true  != val_false);
    REQUIRE(bool_false != val_true );
  }

  SECTION("Check AND operator")
  {
    bool b = val_true  && val_true;
    REQUIRE(b);
    b = val_true && val_false;
    REQUIRE_FALSE(b);
    b = val_false && val_true;
    REQUIRE_FALSE(b);
    b = val_false && val_false;
    REQUIRE_FALSE(b);
  }

  SECTION("Check OR operator")
  {
    bool b = val_true || val_true;
    REQUIRE(b);
    b = val_true || val_false;
    REQUIRE(b);
    b = val_false || val_true;
    REQUIRE(b);
    b = val_false || val_false;
    REQUIRE_FALSE(b);
  }
}


TEST_CASE("test_logical")
{
  TEST_NUMERIC_TYPES(test_logical);
}
