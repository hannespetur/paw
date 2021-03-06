#include "../include/catch.hpp"

#include <deque> // std::deque
#include <list> // std::list
#include <vector> // std::vector

#include <paw/station/algorithm.hpp> // paw::count


using namespace paw;


/*****************************
 * Counting empty containers *
 *****************************/
template <typename T>
void
check_empty_ints()
{
  T ints;
  REQUIRE(paw::count(ints.begin(), ints.end(), 0) == 0);
}


TEST_CASE("Counting empty containers")
{
  SECTION("Empty vector")
  check_empty_ints<std::vector<int> >();

  SECTION("Empty list")
  check_empty_ints<std::list<int> >();

  SECTION("Empty deque")
  check_empty_ints<std::deque<int> >();
}

/*********************
 * Counting elements *
 *********************/
template <typename T>
void
check_int_count_case1()
{
  // 5 zeros total
  T ints = {0, 2, 4, 6, 0, 78, 8, 1, 2, 4, 6, 0, 7, 8, 554, 345, 23, 53, 64, 0, 65, 13, 4444, 0};
  REQUIRE(paw::count(ints.begin(), ints.end(), 0) == 5);
}


TEST_CASE("Counting elements in a vector case 1")
{
  SECTION("Empty vector")
  check_int_count_case1<std::vector<int> >();

  SECTION("Empty list")
  check_int_count_case1<std::list<int> >();

  SECTION("Empty deque")
  check_int_count_case1<std::deque<int> >();
}


template <typename T>
void
check_int_count_case2()
{
  T ints = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  ints.resize(100000, 0);
  REQUIRE(paw::count(ints.begin(), ints.end(), 0u) == 99990);
}


TEST_CASE("Counting elements in a vector case 2")
{
  SECTION("Empty vector")
  check_int_count_case2<std::vector<int> >();
}
