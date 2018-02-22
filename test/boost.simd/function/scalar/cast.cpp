//==================================================================================================
/*!

  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#include <boost/simd/function/scalar/pack_cast.hpp>
#include <scalar_test.hpp>

STF_CASE_TPL( "Check pack_cast return type",  STF_NUMERIC_TYPES )
{
  namespace bs = boost::simd;
  using bs::pack_cast;
  STF_EXPR_IS(pack_cast<std::int8_t  >(T()), std::int8_t );
  STF_EXPR_IS(pack_cast<std::int16_t >(T()), std::int16_t);
  STF_EXPR_IS(pack_cast<std::int32_t >(T()), std::int32_t);
  STF_EXPR_IS(pack_cast<std::int64_t >(T()), std::int64_t);
  STF_EXPR_IS(pack_cast<std::uint8_t >(T()), std::uint8_t);
  STF_EXPR_IS(pack_cast<std::uint16_t>(T()), std::uint16_t);
  STF_EXPR_IS(pack_cast<std::uint32_t>(T()), std::uint32_t);
  STF_EXPR_IS(pack_cast<std::uint64_t>(T()), std::uint64_t);
  STF_EXPR_IS(pack_cast<float        >(T()), float       );
  STF_EXPR_IS(pack_cast<double       >(T()), double      );
}


STF_CASE_TPL( "Check pack_cast behavior",  STF_NUMERIC_TYPES )
{
  namespace bs = boost::simd;
  using bs::pack_cast;
  T x =  2;
  STF_EQUAL(pack_cast<std::int8_t  >(x), static_cast<std::int8_t  >(x));
  STF_EQUAL(pack_cast<std::int16_t >(x), static_cast<std::int16_t >(x));
  STF_EQUAL(pack_cast<std::int32_t >(x), static_cast<std::int32_t >(x));
  STF_EQUAL(pack_cast<std::int64_t >(x), static_cast<std::int64_t >(x));
  STF_EQUAL(pack_cast<std::uint8_t >(x), static_cast<std::uint8_t >(x));
  STF_EQUAL(pack_cast<std::uint16_t>(x), static_cast<std::uint16_t>(x));
  STF_EQUAL(pack_cast<std::uint32_t>(x), static_cast<std::uint32_t>(x));
  STF_EQUAL(pack_cast<std::uint64_t>(x), static_cast<std::uint64_t>(x));
  STF_EQUAL(pack_cast<float        >(x), static_cast<float        >(x));
  STF_EQUAL(pack_cast<double       >(x), static_cast<double       >(x));
}
