#pragma once

#include <iterator> // std::next, std::make_move_iterator
#include <memory> // std::shared_ptr
#include <vector> // std::vector

#include <paw/station/station_options.hpp>


namespace paw
{


template <typename TContainer>
inline
std::vector<std::shared_ptr<TContainer> >
split(TContainer & container, std::size_t const PARTS)
{
  std::vector<std::shared_ptr<TContainer> > split_container;
  split_container.resize(PARTS);
  std::size_t const container_original_size = container.size();

  for (long i = static_cast<long>(PARTS) - 1; i >= 0; --i)
  {
    std::size_t const part_size = container_original_size / PARTS +
                                  (container_original_size % PARTS > static_cast<std::size_t>(i));
    split_container[i] =
      std::make_shared<TContainer>(std::make_move_iterator(std::next(container.end(), -part_size)),
                                   std::make_move_iterator(container.end())
                                   );

    // Make sure the container is smaller now
    if (container.size() > i * (container_original_size / PARTS + 1))
      container.resize(i * (container_original_size / PARTS + 1));
  }

  return split_container;
}


template <typename BidirectionalIterator>
inline
std::vector<std::shared_ptr<std::vector<typename BidirectionalIterator::value_type> > >
split(BidirectionalIterator first, BidirectionalIterator last, StationOptions const & options)
{
  using TContainer = std::vector<typename BidirectionalIterator::value_type>;
  std::size_t const PARTS = options.num_threads;
  std::vector<std::shared_ptr<TContainer> > split_container;
  split_container.resize(PARTS);
  std::size_t const container_original_size = std::distance(first, last);

  for (long i = static_cast<long>(PARTS) - 1; i >= 0; --i)
  {
    std::size_t const part_size = container_original_size / PARTS +
                                  (container_original_size % PARTS > static_cast<std::size_t>(i));
    split_container[i] =
      std::make_shared<TContainer>(std::make_move_iterator(std::next(last, -part_size)),
                                   std::make_move_iterator(last)
                                   );

    last = std::next(last, -part_size);
  }

  return split_container;
}


template <typename BidirectionalIterator>
inline
std::vector<std::shared_ptr<std::vector<typename BidirectionalIterator::value_type> > >
split(BidirectionalIterator first, BidirectionalIterator last, std::size_t const PARTS)
{
  StationOptions options;
  options.set_num_threads(PARTS);
  return split(first, last, options);
}


} // namespace paw
