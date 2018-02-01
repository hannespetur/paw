#pragma once

#include <atomic> // std::atomic
#include <list> // std::list
#include <functional> // std::function


namespace paw
{

class WorkerQueue
{
private:
  /** List of functions in queue. */
  std::list<std::function<void()> > function_queue;

  /** Iterator that points to the current item running.*/
  std::list<std::function<void()> >::iterator item_to_run;

  /** Set when all works/functions have been run.*/
  bool finished = false;

  /** Number of works/functions in queue. */
  std::atomic<std::size_t> queue_size;

  /** Number of completed items run by this queue.*/
  std::size_t completed_items = 0;

public:
  /******************
   * PUBLIC MEMBERS *
   ******************/
  /** Construct a new WorkerQueue instance.*/
  WorkerQueue();

  /** Adds a work/function to the queue.*/
  void add_work_to_queue(std::function<void()> work);

  /** Returns the number of items in queue.*/
  std::size_t get_number_of_items_in_queue() const;

  /** Returns the number of completed items in the queue.*/
  std::size_t get_number_of_completed_items() const;

  /** Runs the queue.*/
  void operator()();

  /** Set the queue to finished.*/
  void finish();

};

} // namespace paw


#ifdef IMPLEMENT_PAW
/* IMPLEMENTATION*/

#include <thread> // std::this_thread::sleep_for

namespace paw
{


WorkerQueue::WorkerQueue()
{
  queue_size = 0;
}


void
WorkerQueue::add_work_to_queue(std::function<void()> work)
{
  // 'function_queue' is std::list because appending to std::list does not invalidate previous
  // iterators.
  function_queue.push_back(work);

  if (queue_size == 0 && completed_items == 0)
    item_to_run = function_queue.begin();

  ++queue_size;
}


std::size_t
WorkerQueue::get_number_of_items_in_queue() const
{
  return queue_size;
}


std::size_t
WorkerQueue::get_number_of_completed_items() const
{
  return completed_items;
}


void
WorkerQueue::operator()()
{
  while (true)
  {
    if (queue_size > 0)
    {
      if (completed_items != 0)
        ++item_to_run;

      (*item_to_run)();
      ++completed_items;
      --queue_size;
    }
    else if (finished)
    {
      return;
    }
    else
    {
      std::this_thread::sleep_for(std::chrono::microseconds(5)); // 0.005 ms
    }
  }
}


void
WorkerQueue::finish()
{
  finished = true;
}


} // namespace paw

#endif //IMPLEMENT_PAW
