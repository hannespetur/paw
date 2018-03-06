#include <numeric>
#include <vector>

#include <boost/simd/algorithm/copy.hpp>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_copy()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(static_cast<T>(2*N+1));
  std::vector<T> ref(static_cast<T>(2*N+3));
  std::vector<T> out(static_cast<T>(2*N+3));
  std::iota(values.begin(), values.end(),T(1));
  std::copy(values.begin(), values.end(), ref.begin());

  // verify we stopped where we should
  REQUIRE(boost::simd::copy(values.data(), values.data() + values.size(), out.data()) == out.data() + values.size());

  // verify values
  REQUIRE(out == ref);
}

TEST_CASE("test copy")
{
  TEST_NUMERIC_TYPES(test_copy);
}
