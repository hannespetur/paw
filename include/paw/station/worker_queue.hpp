#pragma once

#include <atomic> // std::atomic
#include <functional> // std::function
#include <iostream> // std::cout, std::endl
#include <thread> // std::this_thread::sleep_for
#include <list> // std::list


namespace paw
{


class WorkerQueue
{
public:
  // Use std::list because pushing back to lists does not invalidate previous iterators
  std::list<std::function<void()> > function_queue;
  std::list<std::function<void()> >::iterator item_to_run;
  bool finished = false;
  std::atomic<std::size_t> queue_size;
  std::size_t completed_items = 0;


  WorkerQueue();
  void add_work_to_queue(std::function<void()> work);
  std::size_t get_number_of_items_in_queue() const;
  std::size_t get_number_of_completed_items() const;
  void operator()();

};

} // namespace paw


#ifdef IMPLEMENT_PAW
/* IMPLEMENTATION*/

namespace paw
{


WorkerQueue::WorkerQueue()
{
  queue_size = 0;
}


void
WorkerQueue::add_work_to_queue(std::function<void()> work)
{
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


} // namespace paw

#endif //IMPLEMENT_PAW
