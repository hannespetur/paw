#pragma once

#include <algorithm> // std::all_of, ...
#include <atomic> // std::atomic
#include <iterator> // std::iterator_traits

#include <paw/station/internal/algorithm_help_functions.hpp>

#include <paw/station/partition_iterator.hpp> // paw::get_partition_iterators
#include <paw/station/station.hpp> // paw::Station
#include <paw/station/station_options.hpp> // paw::StationOptions

namespace paw
{

/**
 * Checks if unary predicate p returns true for all elements in the range [first, last).
 * \param[in] options Options of the paw::Station.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
bool inline
all_of(StationOptions && options, InputIt first, InputIt last, UnaryPredicate p);


/**
 * Checks if unary predicate p returns true for all elements in the range [first, last).
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
bool inline
all_of(InputIt first, InputIt last, UnaryPredicate p);


/**
 * Checks if unary predicate p returns true for any elements in the range [first, last).
 * \param[in] options Options of the paw::Station.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
bool inline
any_of(StationOptions && options, InputIt first, InputIt last, UnaryPredicate p);


/**
 * Checks if unary predicate p returns true for any elements in the range [first, last).
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
bool inline
any_of(InputIt first, InputIt last, UnaryPredicate p);


/**
 * Returns the number of elements in the range [first, last) satisfying specific criteria.
 * \param[in] options Options of the paw::Station.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] value Value to count.
 */
template <typename InputIt, typename T>
T inline
count(StationOptions && options, InputIt first, InputIt last, T const & value);


/**
 * Returns the number of elements in the range [first, last) satisfying specific criteria.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] value Value to count.
 */
template <typename InputIt, typename T>
T inline
count(InputIt first, InputIt last, T const & value);


/**
 * Returns the number of elements in the range [first, last) satisfying specific criteria.
 * \param[in] options Options of the paw::Station.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
typename std::iterator_traits<InputIt>::difference_type inline
count_if(StationOptions && options, InputIt first, InputIt last, UnaryPredicate p);


/**
 * Returns the number of elements in the range [first, last) satisfying specific criteria.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
typename std::iterator_traits<InputIt>::difference_type inline
count_if(InputIt first, InputIt last, UnaryPredicate p);


/**
 * Assigns the given value to the elements in the range [first, last).
 * \param[in] options Options of the paw::Station.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] value Value to fill.
 */
template <typename InputIt, typename T>
void inline
fill(StationOptions && options, InputIt first, InputIt last, T const & value);


/**
 * Assigns the given value to the elements in the range [first, last).
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] value Value to fill.
 */
template <typename InputIt, typename T>
void inline
fill(InputIt first, InputIt last, T const & value);


/**
 * Applies the given function object f to the result of dereferencing every iterator in the range
 *  [first, last), in order.
 * \param[in] options Options of the paw::Station.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryFunction>
UnaryFunction inline
for_each(StationOptions && options, InputIt first, InputIt last, UnaryFunction p);


/**
 * Applies the given function object f to the result of dereferencing every iterator in the range
 *  [first, last), in order.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryFunction>
UnaryFunction inline
for_each(InputIt first, InputIt last, UnaryFunction p);


/** Checks if unary predicate p returns true for no elements in the range [first, last).
 * \param[in] options Options of the paw::Station.
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
bool inline
none_of(StationOptions && options, InputIt first, InputIt last, UnaryPredicate p);


/** Checks if unary predicate p returns true for no elements in the range [first, last).
 * \param[in] first Beginning of the range.
 * \param[in] last Exclusive end of the range.
 * \param[in] p Predicate to apply to the elements in the range.
 */
template <typename InputIt, typename UnaryPredicate>
bool inline
none_of(InputIt first, InputIt last, UnaryPredicate p);

/**
 * Sorts the elements in the range [first, last) in ascending order.
 * \details The order of equal elements is not guaranteed to be preserved.
 * \param[in] options Options of the paw::Station.
 * \param[in,out] first Beginning of the range.
 * \param[in,out] last Exclusive end of the range.
 */
template <typename InputIt>
void inline
sort(StationOptions && options, InputIt first, InputIt last);


/**
 * Sorts the elements in the range [first, last) in ascending order.
 * \details The order of equal elements is not guaranteed to be preserved.
 * \param[in,out] first Beginning of the range.
 * \param[in,out] last Exclusive end of the range.
 */
template <typename InputIt>
void inline
sort(InputIt first, InputIt last);


} // namespace paw


namespace paw
{

template <typename InputIt, typename UnaryPredicate>
bool inline
all_of(StationOptions && options, InputIt first, InputIt last, UnaryPredicate f)
{
  // Should to be atomic for thread safety
  std::atomic<bool> false_found {
    false
  };

  Station all_of_station(options);
  std::vector<InputIt> partition_iterators = get_partition_iterators(first, last, options);

  for (long i = 0; i < static_cast<long>(partition_iterators.size()) - 1; ++i)
  {
    // If some expression have found to be false we can safely skip the rest of the work
    if (false_found)
      break;

    all_of_station.add_work([&false_found, f](InputIt first, InputIt last)
      {
        if (!std::all_of(first, last, f))
          false_found = true;
      },                       /*function*/
                            partition_iterators[i], /*first*/
                            partition_iterators[i + 1] /*last*/
                            );
  }

  all_of_station.join();
  return !false_found;
}


template <typename InputIt, typename UnaryPredicate>
bool inline
all_of(InputIt first, InputIt last, UnaryPredicate f)
{
  StationOptions options;
  options.chunk_size = 1;
  return paw::all_of(std::move(options), first, last, f);
}


template <typename InputIt, typename UnaryPredicate>
bool inline
any_of(StationOptions && options, InputIt first, InputIt last, UnaryPredicate f)
{
  std::atomic<bool> true_found {
    false
  };                                   // Needs to be atomic for thread safety

  std::vector<InputIt> partition_iterators = get_partition_iterators(first, last, options);
  Station any_of_station(options);

  for (long i = 0; i < static_cast<long>(partition_iterators.size()) - 1; ++i)
  {
    if (true_found)
      break;

    any_of_station.add_work([&true_found, f](InputIt first, InputIt last)
      {
        if (std::any_of(first, last, f))
          true_found = true;
      } /*function*/,
                            partition_iterators[i], /*first*/
                            partition_iterators[i + 1] /*last*/
                            );
  }

  any_of_station.join();
  return true_found;
}


template <typename InputIt, typename UnaryPredicate>
bool inline
any_of(InputIt first, InputIt last, UnaryPredicate f)
{
  StationOptions options;
  options.chunk_size = 1;
  return paw::any_of(std::move(options), first, last, f);
}


template <typename InputIt, typename T>
T inline
count(StationOptions && options, InputIt first, InputIt last, T const & value)
{
  std::vector<InputIt> partition_iterators = get_partition_iterators(first, last, options);
  std::vector<std::shared_ptr<T> > counts;
  Station count_station(options);

  for (long i = 0; i < static_cast<long>(partition_iterators.size()) - 1; ++i)
  {
    counts.push_back(std::make_shared<T>(0));
    count_station.add_work([value](InputIt first, InputIt last, std::shared_ptr<T> ret)
      {
        *ret = std::count(first, last, value);
      } /*function*/,
                           partition_iterators[i], /*first*/
                           partition_iterators[i + 1], /*last*/
                           counts.back()
                           );
  }

  // Merge region
  T sum = 0;
  count_station.join(); // Join all threads

  for (auto const & count : counts)
    sum += *count;

  return sum;
}


template <typename InputIt, typename T>
T inline
count(InputIt first, InputIt last, T const & value)
{
  StationOptions options;
  options.chunk_size = 0; // Partition the container evenly
  return count(std::move(options), first, last, value);
}


template <typename InputIt, typename UnaryPredicate>
typename std::iterator_traits<InputIt>::difference_type inline
count_if(StationOptions && options, InputIt first, InputIt last, UnaryPredicate p)
{
  using T = typename std::iterator_traits<InputIt>::difference_type;
  std::vector<InputIt> partition_iterators = get_partition_iterators(first, last, options);
  std::vector<std::shared_ptr<T> > counts;
  Station count_if_station(options);

  for (long i = 0; i < static_cast<long>(partition_iterators.size()) - 1; ++i)
  {
    counts.push_back(std::make_shared<T>(0));
    count_if_station.add_work([p](InputIt first, InputIt last, std::shared_ptr<T> ret)
      {
        *ret = std::count_if(first, last, p);
      } /*function*/,
                              partition_iterators[i], /*first*/
                              partition_iterators[i + 1], /*last*/
                              counts.back()
                              );
  }

  // Merge region
  T sum = 0;
  count_if_station.join();

  for (auto const & count : counts)
    sum += *count;

  return sum;
}


template <typename InputIt, typename UnaryPredicate>
typename std::iterator_traits<InputIt>::difference_type inline
count_if(InputIt first, InputIt last, UnaryPredicate p)
{
  return paw::count_if(StationOptions(), first, last, p);
}


template <typename InputIt, typename T>
void inline
fill(StationOptions && options, InputIt first, InputIt last, T const & value)
{
  std::vector<InputIt> partition_iterators = get_partition_iterators(first, last, options);
  std::vector<std::shared_ptr<T> > counts;
  Station fill_station(options);

  for (long i = 0; i < static_cast<long>(partition_iterators.size()) - 1; ++i)
  {
    fill_station.add_work([value](InputIt first, InputIt last)
      {
        std::fill(first, last, value);
      } /*function*/,
                          partition_iterators[i], /*first*/
                          partition_iterators[i + 1] /*last*/
                          );
  }

  fill_station.join();
}


template <typename InputIt, typename T>
void inline
fill(InputIt first, InputIt last, T const & value)
{
  StationOptions options;
  options.chunk_size = 0; // Partition evenly
  fill(std::move(options), first, last, value);
}


template <typename InputIt, typename UnaryFunction>
UnaryFunction inline
for_each(StationOptions && options, InputIt first, InputIt last, UnaryFunction f)
{
  std::vector<InputIt> partition_iterators = get_partition_iterators(first, last, options);
  Station for_each_station(options);

  for (long i = 0; i < static_cast<long>(partition_iterators.size()) - 1; ++i)
    for_each_station.add_work([f](InputIt it){
        f(*it);
      } /*function*/, partition_iterators[i] /*it*/);

  for_each_station.join();
  return f;
}


template <typename InputIt, typename UnaryFunction>
UnaryFunction inline
for_each(InputIt first, InputIt last, UnaryFunction f)
{
  StationOptions options;
  options.chunk_size = 1;
  return paw::for_each(std::move(options), first, last, f);
}


template <typename InputIt, typename UnaryPredicate>
bool inline
none_of(StationOptions && options, InputIt first, InputIt last, UnaryPredicate f)
{
  return !paw::any_of(options, first, last, f);
}


template <typename InputIt, typename UnaryPredicate>
bool inline
none_of(InputIt first, InputIt last, UnaryPredicate f)
{
  return !paw::any_of(first, last, f);
}


template <typename InputIt>
void inline
sort(StationOptions && options, InputIt first, InputIt last)
{
  std::vector<InputIt> partition_iterators = get_partition_iterators(first, last, options);
  Station sort_station(options);

  for (long i = 0; i < static_cast<long>(partition_iterators.size()) - 1; ++i)
  {
    auto work = [](InputIt first, InputIt last){
                  std::sort(first, last);
                };
    sort_station.add_work(std::move(work),
                          partition_iterators[i], /*first*/
                          partition_iterators[i + 1] /*last*/
                          );
  }

  // After this join, all partitions are sorted, then we need to merge the sorted partition
  sort_station.join();

  // I have tried to implement a multi-threaded merge, but actually the single threaded version
  // was always faster! x(
  for (std::size_t d = 2; d < partition_iterators.size(); ++d)
  {
    std::inplace_merge(partition_iterators[0], partition_iterators[d - 1], partition_iterators[d]);
  }
}


template <typename InputIt>
void inline
sort(InputIt first, InputIt last)
{
  StationOptions options;
  std::size_t const container_size = std::distance(first, last);

  if (options.num_threads > 2 && container_size >= 1000 && container_size <= 10000)
  {
    options.set_num_threads(2);
  }
  else if (options.num_threads > 4 && container_size >= 10000 && container_size <= 100000)
  {
    options.set_num_threads(4);
  }

  paw::sort(std::move(options), first, last);
}


} // namespace paw
