
//             Copyright NumScale SAS 2017.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <numeric>
#include <vector>

#include <boost/simd/algorithm/equal.hpp>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"

using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_equal()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values1(2 * N + 1), values2(2 * N + 1);;
  std::iota(values1.begin(), values1.end(), T(1));
  std::iota(values2.begin(), values2.end(), T(1));

   REQUIRE(std::equal(values1.begin(), values1.end(), values2.begin()) ==
           boost::simd::equal(values1.data(), values1.data()+values2.size(), values2.data())
           );

   values1[N] = T(0);

   REQUIRE(std::equal(values1.begin(), values1.end(), values2.begin()) ==
           boost::simd::equal(values1.data(), values1.data()+values1.size(), values2.data())
           );
}

TEST_CASE("test equal")
{
  TEST_NUMERIC_TYPES(test_equal);
}
