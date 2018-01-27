#include <algorithm>
#include <bitset>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include <paw/station/internal/data_simulation.hpp>

#include <paw/station/algorithm.hpp>
#include <paw/station/join.hpp>
#include <paw/station/split.hpp>
#include <paw/station/station.hpp>
#include <paw/station/worker_queue.hpp>


int
main(int argc, char ** argv)
{
  srand(42);  // Seed is 42

  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <NUM_INTS>" << std::endl;
    std::exit(1);
  }

  int const num_ints = std::stoi(argv[1]);
  std::vector<int> ints = paw::stations_internal::get_random_ints<std::vector<int> >(num_ints);
  auto start = std::chrono::system_clock::now();
  paw::sort(ints.begin(), ints.end());
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Total number of ints sorted are " << ints.size() << std::endl;
  std::cout << "Sorting duration was " << diff.count() << " seconds." << std::endl;

  // Make sure the vector is sorted
  for (auto it = ints.begin() + 1; it != ints.end(); ++it)
  {
    if (*it < *(it - 1))
      std::cout << "NOT SORTED" << std::endl;
  }
}
