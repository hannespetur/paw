#include <algorithm>
#include <numeric>
#include <vector>

#include <boost/simd/algorithm/count_if.hpp>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_count_if()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(2 * N + 1);
  std::iota(values.begin(), values.end(),T(1));

  auto c = std::count_if(values.begin(), values.end(), is_odd);
  auto bc = boost::simd::count_if(values.data(), values.data()+2*N+1, is_odd);

  REQUIRE(bc == c);
}

TEST_CASE("test count_if")
{
  TEST_NUMERIC_TYPES(test_count_if);
}
