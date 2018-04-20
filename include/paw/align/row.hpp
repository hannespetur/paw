#pragma once

#include <cstdint>
#include <vector>

#include <boost/simd/pack.hpp>
#include <boost/simd/memory/allocator.hpp>


namespace paw
{


template <typename Tuint>
struct Row1
{
  using uint_t = Tuint;
  using Tpack = boost::simd::pack<Tuint>;
  using Tlogical_pack = boost::simd::pack<boost::simd::logical<Tuint> >;
  using Vectors = std::vector<Tpack, boost::simd::allocator<Tpack> >;

  std::size_t static constexpr vector_size = boost::simd::cardinal_of<Tpack>();

  std::size_t const n_elements;
  Vectors vectors;

  /* CONSTRUCTORS */
  Row1(std::size_t const _n_elements);
  Row1(std::size_t const _n_elements, Tuint const val);

};


} // namespace paw


#if defined(IMPLEMENT_PAW)

#include <iomanip> // std::setw


namespace paw
{

template<typename Tuint>
std::size_t constexpr Row1<Tuint>::vector_size;


template <typename Tuint>
Row1<Tuint>::Row1(std::size_t const _n_elements)
  : n_elements(_n_elements)
{
  Tpack my_vector;
  vectors.resize((n_elements + vector_size - 1) / vector_size, my_vector);
}


template <typename Tuint>
Row1<Tuint>::Row1(std::size_t const _n_elements, Tuint const val)
  : n_elements(_n_elements)
{
  Tpack my_vector {
    val
  };

  vectors.resize((n_elements + vector_size - 1) / vector_size, my_vector);
}


template <typename Tint>
std::ostream &
operator<<(std::ostream & ss, Row1<Tint> const & r)
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
