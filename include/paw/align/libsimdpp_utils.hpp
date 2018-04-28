#pragma once

#include <array>
#include <cassert>
#include <vector>

#include <simdpp/simd.h>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

struct Row;

#if SIMDPP_USE_AVX2
constexpr int S = 32;
#elif SIMDPP_USE_SSE2
constexpr int S = 16;
#else
constexpr int S = 16;
#endif


namespace T
{

#define PAW_USE_UINT8

#if defined(PAW_USE_UINT8)
using uint = uint8_t;
using pack = simdpp::uint8<S / sizeof(uint), void>;
using mask = simdpp::mask_int8<S / sizeof(uint), void>;
#else
using uint = uint16_t;
using pack = simdpp::uint16<S / sizeof(uint), void>;
using mask = simdpp::mask_int16<S / sizeof(uint), void>;
#endif

using vec_pack = std::vector<pack, simdpp::aligned_allocator<pack, sizeof(pack)> >;
using vec_uint = std::vector<uint, simdpp::aligned_allocator<uint, sizeof(uint)> >;
using arr_row = std::array<Row, 4>;
using arr_uint = std::array<uint, S / sizeof(uint)>;

} // namespace T


struct Row
{
  long const n_elements = 0;
  T::vec_pack vectors;

  /* CONSTRUCTORS */
  Row(std::size_t const _n_elements)
    : n_elements(_n_elements)
    , vectors(0)
  {
    T::pack my_vector = simdpp::make_zero();
    vectors.resize((n_elements + T::pack::length - 1) / T::pack::length, my_vector);
  }


  Row(std::size_t const _n_elements, T::uint const val)
    : n_elements(_n_elements)
  {
    T::pack my_vector = simdpp::make_uint(val);
    vectors.resize((n_elements + T::pack::length - 1) / T::pack::length, my_vector);
  }
};


inline T::pack
shift_one_right(T::pack pack, T::uint const left)
{
  std::array<T::uint, T::pack::length + 1> vec;
  //vec.fill(left);
  vec[0] = left;
  simdpp::store_u(&vec[1], pack);
  return simdpp::load_u(&vec[0]);
}


inline T::pack
shift_one_right(T::pack pack, T::uint const left, std::array<long, S / sizeof(T::uint)> const & reductions)
{
  std::array<T::uint, T::pack::length + 1> vec;
  vec[0] = left;
  simdpp::store_u(&vec[1], pack);

  for (long e = 1; e < static_cast<long>(T::pack::length); ++e)
  {
    long const val = static_cast<long>(vec[e]) + reductions[e - 1] - reductions[e];
    vec[e] = val > 0 ? val : 0;
  }

  return simdpp::load_u(&vec[0]);
  //#elif SIMDPP_USE_SSE2
  //return simdpp::align16<15, 16>(static_cast<T::pack>(simdpp::make_uint(left)), pack);
  //#else
  //return simdpp::align16<15, 16>(static_cast<T::pack>(simdpp::make_uint(left)), pack);
  //#endif
}


inline T::pack
shift_one_right(T::pack pack)
{
  return shift_one_right(pack, 0);
}


template<typename Tpack>
void
print_pack(Tpack const & pack)
{
  using T = typename Tpack::uint_element_type;

  // Guard for when the pack is empty
  if (pack.length == 0)
  {
    std::cout << "()";
    return;
  }

  std::vector<T, simdpp::aligned_allocator<T, sizeof(T)> > vector;
  vector.resize(pack.length);
  simdpp::store_u(&vector[0], pack);

  std::cout << "(" << static_cast<long>(vector[0]);

  for (long i = 1; i < static_cast<long>(vector.size()); ++i)
  {
    std::cout << "," << static_cast<long>(vector[i]);
  }

  std::cout << ")" << std::endl;
}


inline void
print_score_vector_standard(Row const & vX)
{
  long const t = vX.vectors.size();

  if (t == 0)
    return;

  T::vec_uint vec(vX.vectors[0].length, 0);
  std::vector<T::vec_uint> m(vX.vectors.size(), vec);

  for (long v = 0; v < t; ++v)
  {
    simdpp::store_u(&m[v][0], vX.vectors[v]);
  }

  for (long j = 0; j < vX.n_elements; ++j)
  {
    std::size_t const v = j % t;
    std::size_t const e = j / t;
    assert(v < m.size());
    assert(e < m[v].size());
    std::cout << std::setw(4) << static_cast<uint64_t>(m[v][e]) << " ";
  }

  std::cout << "(" << t << " vectors, " << vX.n_elements << " elements)\n";
}



inline void
print_score_vectors(Row const & vH,
                    Row const & vH_up,
                    Row const & vE,
                    Row const & vF,
                    Row const & vF_up,
                    Row const & vW
                    )
{
  std::cout << "Standard H_up  : "; print_score_vector_standard(vH_up);
  std::cout << "Standard H     : "; print_score_vector_standard(vH);
  std::cout << "Standard E     : "; print_score_vector_standard(vE);
  std::cout << "Standard F_up  : "; print_score_vector_standard(vF_up);
  std::cout << "Standard F     : "; print_score_vector_standard(vF);
  std::cout << "Standard W     : "; print_score_vector_standard(vW);
  std::cout << "=====\n";
}


} //namespace SIMDPP_ARCH_NAMESPACE
} // anon namespace
