#pragma once

#include <memory> // std::unique_ptr
#include <thread> // std::thread
#include <utility> // std::forward

#include <paw/station/station_options.hpp> // paw::station::WorkerQueue
#include <paw/station/partition_iterator.hpp> // paw::station::get_partition_iterators
#include <paw/station/worker_queue.hpp> // paw::station::WorkerQueue


/****************
 * DECLERATIONS *
 ****************/

/** \file paw/station/station.hpp
 * Main header file for paw Station.
 */

/** Top-level namespace of the paw library.*/
namespace paw
{

class Station
{
public:
  using Queues = std::vector<std::unique_ptr<WorkerQueue> >;

private:
  /*********************
  * STATION STATISTICS *
  *********************/
  /** Set if the station has joined its worker queues, i.e. all threads have finished running. */
  bool joined = false;

  /** Number of jobs the main thread has run. */
  std::size_t main_thread_work_count = 0;

  /** std::thread instances in this station. */
  std::vector<std::thread> workers;

  /** Each thread has a unique worker queue. */
  Queues queues;

public:
  /** Options and policies this station will follow */
  StationOptions options;

  /** Constructs a new Station instance.
   * \param[in] _options Station options to use.
   */
  Station(StationOptions _options);

  /** Constructs a new Station instance.
   * \param[in] num_threads Number of threads to use (includes the main thread).
   * \param[in] max_queue_size Maximum number of elements in each queue.
   */
  Station(std::size_t const num_threads = 1, std::size_t const max_queue_size = 2);

  /** Deconstructs the Station. */
  ~Station();

  /***************************
  * GENERAL MEMBER FUNCTIONS *
  ***************************/
  template <typename TWork, typename ... Args>
  void inline
  add_work(TWork && work, Args ... args);

  template <typename TWork, typename ... Args>
  void inline
  add_work_with_thread_id(TWork && work, Args ... args);

  template <typename TWork, typename ... Args>
  void inline
  add(TWork && work, Args ... args);

  template <typename TWork, typename ... Args>
  void inline
  add_to_thread(std::size_t const thread_id, TWork && work, Args ... args);

  template <typename TWork, typename ... Args>
  void inline
  add_to_thread_with_thread_id(std::size_t const thread_id, TWork && work, Args ... args);

  Queues::const_iterator find_smallest_queue(std::size_t & smallest_size);

  std::string join();

  /***************************
  * PRIVATE MEMBER FUNCTIONS *
  ****************************/
private:
  void resize_queues_and_workers(std::size_t const new_size);

};


/*************
 * TEMPLATES *
 *************/
template <typename TWork, typename ... Args>
void inline
Station::add_work(TWork && work, Args ... args)
{
  if (workers.size() == 0)
  {
    work(args ...);
    ++main_thread_work_count;
  }
  else
  {
    // We have some workers, so let's assign the work to the smallest worker queue
    std::size_t smallest_size = -1;

    // If we have any worker threads, check who has the smallest queue
    auto min_queue_it = find_smallest_queue(smallest_size);

    if (smallest_size < options.max_queue_size)
    {
      (*min_queue_it)->add_work_to_queue(std::bind(std::forward<TWork>(work), args ...));
    }
    else
    {
      work(args ...);   // If all queues are of maximum size, use the boss thread instead
      ++main_thread_work_count;
    }
  }
}


template <typename TWork, typename ... Args>
void inline
Station::add_work_with_thread_id(TWork && work, Args ... args)
{
  if (workers.size() == 0)
  {
    work(options.num_threads - 1, args ...);
    ++main_thread_work_count;
  }
  else
  {
    // We have some workers, so let's assign the work to the smallest worker queue
    std::size_t smallest_size = -1;

    // If we have any worker threads, check who has the smallest queue
    Queues::const_iterator min_queue_it = find_smallest_queue(smallest_size);

    if (smallest_size < options.max_queue_size)
    {
      long const thread_id = std::distance(queues.cbegin(), min_queue_it);
      (*min_queue_it)->add_work_to_queue(std::bind(std::forward<TWork>(work), thread_id, args ...));
    }
    else
    {
      std::size_t const thread_count = options.num_threads;
      work(thread_count - 1, args ...);   // If all queues are of maximum size, use the boss thread instead
      ++main_thread_work_count;
    }
  }
}


template <typename TWork, typename ... Args>
void inline
Station::add(TWork && work, Args ... args)
{
  // For backwards compability
  this->add_work(std::forward<TWork>(work), args ...);
}


template <typename TWork, typename ... Args>
void inline
Station::add_to_thread(std::size_t const thread_id, TWork && work, Args ... args)
{
  std::size_t const thread_count = options.num_threads;

  if (thread_id % thread_count == thread_count - 1)
  {
    work(args ...);
    ++main_thread_work_count;
  }
  else
  {
    queues[thread_id % thread_count]->add_work_to_queue(std::bind(std::forward<TWork>(work),
                                                                  args ...)
                                                        );
  }
}


} // namespace paw


#if defined(IMPLEMENT_PAW) || defined(__JETBRAINS_IDE__)

#include <iostream> // std::cout

namespace paw
{

Station::Station(StationOptions _options)
  : options(_options)
{
  resize_queues_and_workers(options.num_threads - 1);
}


Station::Station(std::size_t const num_threads, std::size_t const max_queue_size)
{
  options.num_threads = num_threads;
  options.max_queue_size = max_queue_size;
  resize_queues_and_workers(options.num_threads - 1);
}


Station::~Station()
{
  if (not joined)
    join();
}


Station::Queues::const_iterator
Station::find_smallest_queue(std::size_t & smallest_size)
{
  Queues::const_iterator smallest_q_it;

  for (auto q_it = queues.cbegin(); q_it != queues.cend(); ++q_it)
  {
    std::size_t const current_queue_size = (*q_it)->get_number_of_items_in_queue();

    if (current_queue_size < smallest_size)
    {
      smallest_size = current_queue_size;
      smallest_q_it = q_it;

      // Check if any queue is empty, and if so then we don't need to look any further
      if (smallest_size == 0)
        break;
    }
  }

  return smallest_q_it;
}


std::string
Station::join()
{
  std::ostringstream ss;
  ss << main_thread_work_count;

  for (int i = 0; i < static_cast<int>(options.num_threads) - 1; ++i)
  {
    queues[i]->finish();
    workers[i].join();
    ss << "/" << queues[i]->get_number_of_completed_items();
  }

  joined = true;
  return ss.str();
}


void
Station::resize_queues_and_workers(std::size_t const new_size)
{
  for (std::size_t i = 0; i < new_size; ++i)
  {
    queues.push_back(std::unique_ptr<WorkerQueue>(new WorkerQueue()));
    workers.push_back(std::thread(std::ref(*queues[i])));
  }
}


} // namespace paw

#endif // IMPLEMENT_PAW
