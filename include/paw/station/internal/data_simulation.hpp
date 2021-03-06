#pragma once

#include <cstdlib>
#include <vector> // std::vector

namespace paw
{
namespace stations_internal
{

template <typename TContainer>
TContainer inline
get_random_ints(std::size_t const N)
{
  TContainer ints;

  for (std::size_t i = 0; i < N; ++i)
  {
    // Random number between (and including) -1,000,000 and 1,000,000
    ints.push_back(rand() % 200000001 - 100000000);
  }

  return ints;
}


// std::vector specialization
template <>
std::vector<int> inline
get_random_ints(std::size_t const N)
{
  std::vector<int> ints;
  ints.reserve(N);

  for (std::size_t i = 0; i < N; ++i)
  {
    // Random number between (and including) -1,000,000 and 1,000,000
    ints.push_back(rand() % 200000001 - 100000000);
  }

  return ints;
}


} // namespace stations_internal
} // namespace paw
