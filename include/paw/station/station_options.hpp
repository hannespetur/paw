#pragma once

#include <thread> // std::thread::hardware_concurrency()


namespace paw
{

class StationOptions
{
  friend class Station; /** Allow stations to see your privates. */

public:
  /** Maximum size of the worker queues, including jobs which are running */
  std::size_t max_queue_size = 2;

  /** Number of items in each chunk of work to process. If 0, then the work will be evenly
   * distributed among all threads.
   */
  std::size_t chunk_size = 0;

  /** Number of threads to use, including the main thread */
  std::size_t num_threads = std::thread::hardware_concurrency();

  /** 0 is quite mode, 1 can output warnings, and 2 will output warnings and statistics to
   *  std::cout.
   */
  std::size_t verbosity = 0;

  void set_num_threads(std::size_t const _num_threads);
};

} // namespace paw


#if defined(IMPLEMENT_PAW) || defined(__JETBRAINS_IDE__)

namespace paw
{

void
StationOptions::set_num_threads(std::size_t const _num_threads)
{
  if (_num_threads == 0)
    num_threads = std::max(static_cast<uint32_t>(1), std::thread::hardware_concurrency());
  else
    num_threads = _num_threads;
}


} // namespace paw

#endif // IMPLEMENT_PAW
