#include <numeric> // std::iota

#include <boost/align/aligned_allocator.hpp>
#include <boost/simd/pack.hpp>
#include <boost/simd/function/abs.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/store.hpp>

#include "../test_numeric_types.hpp"
#include "../../../include/catch.hpp"


namespace
{

template<typename T>
T
absolute(T const a)
{
  if (a >= T(0))
    return a;
  else
    return -a;
}


template<typename T>
void
test_abs()
{
  namespace ba = boost::alignment;
  namespace bs = boost::simd;

  static const int N = boost::simd::pack<T>::static_size;
  static const int M = 2 * N + 1;

  using pack_t = boost::simd::pack<T>;
  using vec_t = std::vector<T, ba::aligned_allocator<T, pack_t::alignment> >;

  vec_t out(static_cast<T>(M));
  vec_t ref(M);

  std::iota(out.begin(), out.end(), T(-3));
  std::iota(ref.begin(), ref.end(), T(-3));

  for (int i = 0; i < M; i += N)
  {
    pack_t v0 = bs::load<pack_t>(&out[i]);
    bs::store(bs::abs(v0), &out[i]);
  }

  std::transform(ref.begin(), ref.end(), ref.begin(), absolute<T>);

  REQUIRE(out == ref);
}


TEST_CASE("test abs")
{
  test_abs<uint8_t>();
  test_abs<uint16_t>();
  test_abs<uint32_t>();
  test_abs<uint64_t>();

  test_abs<int8_t>();
  test_abs<int16_t>();
  test_abs<int32_t>();
  test_abs<int64_t>();

  test_abs<float>();
  test_abs<double>();
}

} // anon namespace
