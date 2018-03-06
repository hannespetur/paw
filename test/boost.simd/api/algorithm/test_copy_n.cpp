#include <numeric>
#include <vector>

#include <boost/simd/algorithm/copy_n.hpp>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"


using namespace boost::simd;
using namespace boost::alignment;


template<typename T>
void
test_copy_n()
{
  static const int N = pack<T>::static_size;

  std::vector<T> values(static_cast<T>(2 * N + 1));
  std::vector<T> ref(static_cast<T>(2 * N + 3));
  std::vector<T> out(static_cast<T>(2 * N + 3));
  std::iota(values.begin(), values.end(),T(1));

  std::copy(values.begin(), values.begin() + N, ref.begin());

  // verify we stopped where we should
  auto s = out.data() + N;
  auto bs = boost::simd::copy_n(values.data(), N, out.data());
  REQUIRE(bs == s);

  // verify values
  REQUIRE(out == ref);
}

TEST_CASE("test copy_n")
{
  TEST_NUMERIC_TYPES(test_copy_n);
}
