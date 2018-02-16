#pragma once

#include <cstdint>

#include <boost/simd.hpp>


namespace paw
{


template <typename Tuint>
struct Row
{
  using uint_t = Tuint;
  using Tpack = boost::simd::pack<Tuint>;
  using Vectors = std::vector<Tpack, boost::simd::allocator<Tpack> >;

  std::size_t static constexpr vector_size = boost::simd::cardinal_of<Tpack>();

  std::size_t const n_elements;
  Vectors vectors;

  /* CONSTRUCTORS */
  Row(std::size_t const _n_elements);
  Row(std::size_t const _n_elements, Tuint const val);
};


} // namespace paw


#ifdef IMPLEMENT_PAW

#include <iomanip> // std::setw

namespace paw
{

template <typename Tuint>
Row<Tuint>::Row(std::size_t const _n_elements)
  : n_elements(_n_elements)
{
  Tpack my_vector;
  vectors.resize((n_elements + vector_size - 1) / vector_size, my_vector);
}


template <typename Tuint>
Row<Tuint>::Row(std::size_t const _n_elements, Tuint const val)
  : n_elements(_n_elements)
{
  Tpack my_vector {
    val
  };
  vectors.resize((n_elements + vector_size - 1) / vector_size, my_vector);
}


template <typename Tint>
std::ostream &
operator<<(std::ostream & ss, Row<Tint> const & r)
{
  std::size_t j = 0;

  for (std::size_t e = 0; e < r.vector_size; ++e)
  {
    for (std::size_t v = 0; v < r.vectors.size(); ++v)
    {
      ss << std::setw(6) << static_cast<int64_t>(r.vectors[v][e]);
      ++j;

      if (j == r.n_elements)
        return ss;
    }
  }

  return ss;
}


} // namespace paw


#endif // IMPLEMENT_PAW
