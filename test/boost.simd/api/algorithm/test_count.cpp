#include <algorithm>
#include <numeric>
#include <vector>

#include <boost/simd/algorithm/count.hpp>

#include "../test_numeric_types.hpp"
#include <catch.hpp>


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_count()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(static_cast<T>(2 * N + 1));
  std::iota(values.begin(), values.end(), T(1));

  auto c = std::count(values.begin(), values.end(), T(5));
  auto bc = boost::simd::count(values.data(), values.data()+2*N+1, T(5));

  REQUIRE(bc == c);
}

TEST_CASE("test count")
{
  TEST_NUMERIC_TYPES(test_count);
}
