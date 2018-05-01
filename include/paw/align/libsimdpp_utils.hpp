#pragma once

#include <array>
#include <cassert>
#include <vector>

#include <simdpp/simd.h>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

//struct Row;
struct Row8;
struct Row16;

#if SIMDPP_USE_AVX512BW
constexpr int S = 64;
#elif SIMDPP_USE_AVX2
constexpr int S = 32;
#elif SIMDPP_USE_SSE2
constexpr int S = 16;
#else
constexpr int S = 16;
#endif


namespace T
{

#define PAW_USE_UINT8

using pack_8 = simdpp::uint8<S / sizeof(uint8_t), void>;
using pack_16 = simdpp::uint16<S / sizeof(uint16_t), void>;

#if defined(PAW_USE_UINT8)
using pack = pack_8;
using row = Row8;
#else
using pack = pack_16;
using row = Row16;
#endif

using mask_8 = pack_8::mask_vector_type;
using uint_8 = pack_8::uint_element_type;
using vec_pack_8 = std::vector<pack_8, simdpp::aligned_allocator<pack_8, sizeof(pack_8)> >;
using vec_uint_8 = std::vector<uint_8, simdpp::aligned_allocator<uint_8, sizeof(uint_8)> >;
using arr_row_8 = std::array<Row8, 4>;
using arr_uint_8 = std::array<uint_8, S / sizeof(uint_8)>;

using mask_16 = pack_16::mask_vector_type;
using uint_16 = pack_16::uint_element_type;
using vec_pack_16 = std::vector<pack_16, simdpp::aligned_allocator<pack_16, sizeof(pack_16)> >;
using vec_uint_16 = std::vector<uint_16, simdpp::aligned_allocator<uint_16, sizeof(uint_16)> >;
using arr_row_16 = std::array<Row16, 4>;
using arr_uint_16 = std::array<uint_16, S / sizeof(uint_16)>;

using mask = pack::mask_vector_type;
using uint = pack::uint_element_type;
using vec_pack = std::vector<pack, simdpp::aligned_allocator<pack, sizeof(pack)> >;
using vec_uint = std::vector<uint, simdpp::aligned_allocator<uint, sizeof(uint)> >;
using arr_row = std::array<row, 4>;
using arr_uint = std::array<uint, S / sizeof(uint)>;

} // namespace T


struct Row8
{
  long const n_elements = 0;
  T::vec_pack_8 vectors;

  /* CONSTRUCTORS */
  Row8(std::size_t const _n_elements)
    : n_elements(_n_elements)
    , vectors(0)
  {
    T::pack_8 my_vector = simdpp::make_zero();
    vectors.resize((n_elements + T::pack_8::length - 1) / T::pack_8::length, my_vector);
  }


  Row8(std::size_t const _n_elements, T::uint_8 const val)
    : n_elements(_n_elements)
  {
    T::pack_8 my_vector = simdpp::make_uint(val);
    vectors.resize((n_elements + T::pack_8::length - 1) / T::pack_8::length, my_vector);
  }
};


struct Row16
{
  long const n_elements = 0;
  T::vec_pack_16 vectors;

  /* CONSTRUCTORS */
  Row16(std::size_t const _n_elements)
    : n_elements(_n_elements)
    , vectors(0)
  {
    T::pack_16 my_vector = simdpp::make_zero();
    vectors.resize((n_elements + T::pack_16::length - 1) / T::pack_16::length, my_vector);
  }


  Row16(std::size_t const _n_elements, T::uint_16 const val)
    : n_elements(_n_elements)
  {
    T::pack_16 my_vector = simdpp::make_uint(val);
    vectors.resize((n_elements + T::pack_16::length - 1) / T::pack_16::length, my_vector);
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
print_score_vector_standard(T::row const & vX)
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
print_score_vectors(T::row const & vH,
                    T::row const & vH_up,
                    T::row const & vE,
                    T::row const & vF,
                    T::row const & vF_up,
                    T::row const & vW
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
