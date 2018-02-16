//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_CONSTANT_NBEXPONENTBITS_HPP_INCLUDED
#define BOOST_SIMD_CONSTANT_NBEXPONENTBITS_HPP_INCLUDED

/*!
  @ingroup group-constant
  @defgroup constant-Nbexponentbits Nbexponentbits (function template)

  Generates a constant representing the number of exponent bits of a floating point type.

  @headerref{<boost/simd/constant/nbexponentbits.hpp>}

  @par Description

  1.  @code
      template<typename T> as_integer_t<T> Nbexponentbits();
      @endcode

  2.  @code
      template<typename T> as_integer_t<T> Nbexponentbits( boost::simd::as_<T> const& target );
      @endcode

    Generates a value of type `as_integer_t<T>` that evaluates to the number of bits used to
    represents the exponent of an IEEE754 floating-point value.

  @par Parameters

  | Name                | Description                                                         |
  |--------------------:|:--------------------------------------------------------------------|
  | **target**          | a [placeholder](@ref type-as) value encapsulating the constant type |

  @par Return Value
  A value of type @c as_integer_t<T> that evaluates to:

  | Type        | double      | float         |
  |:------------|:------------|---------------|
  | **Values**  |   11        |      8        |

  @par Requirements
  - **T** models IEEEValue
**/

#include <boost/simd/constant/scalar/nbexponentbits.hpp>
#include <boost/simd/constant/simd/nbexponentbits.hpp>

#endif
