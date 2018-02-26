#include <numeric>
#include <vector>

#include <boost/simd/algorithm.hpp>
#include <boost/align/aligned_allocator.hpp>

#include "../../include/catch.hpp"

#include "algorithm/test_all_of.hpp"
#include "algorithm/test_any_of.hpp"
#include "algorithm/test_copy.hpp"
#include "algorithm/test_copy_n.hpp"
#include "algorithm/test_count.hpp"
#include "algorithm/test_count_if.hpp"
#include "algorithm/test_equal.hpp"
#include "algorithm/test_fill.hpp"
#include "algorithm/test_find.hpp"
#include "algorithm/test_find_if.hpp"
#include "algorithm/test_find_if_not.hpp"


template<typename T>
void
test_algorithm()
{
  test_all_of<T>();
  test_any_of<T>();
  test_copy<T>();
  test_copy_n<T>();
  test_count<T>();
  test_count_if<T>();
  test_equal<T>();
  test_fill<T>();
  test_find<T>();
  test_find_if<T>();
  test_find_if_not<T>();
}


TEST_CASE("Check algorithms with 8 bit unsigned integers.")
{
  test_algorithm<uint8_t>();
}


TEST_CASE("Check algorithms with 16 bit unsigned integers.")
{
  test_algorithm<uint16_t>();
}


TEST_CASE("Check algorithms with 32 bit unsigned integers.")
{
  test_algorithm<uint32_t>();
}


TEST_CASE("Check algorithms with 64 bit unsigned integers.")
{
  test_algorithm<uint64_t>();
}


TEST_CASE("Check algorithms with 8 bit integers.")
{
  test_algorithm<int8_t>();
}


TEST_CASE("Check algorithms with 16 bit integers.")
{
  test_algorithm<int16_t>();
}


TEST_CASE("Check algorithms with 32 bit integers.")
{
  test_algorithm<int32_t>();
}


TEST_CASE("Check algorithms with 64 bit integers.")
{
  test_algorithm<int64_t>();
}


TEST_CASE("Check algorithms with single-precision numbers.")
{
  test_algorithm<float>();
}


TEST_CASE("Check algorithms with double-precision numbers.")
{
  test_algorithm<double>();
}
