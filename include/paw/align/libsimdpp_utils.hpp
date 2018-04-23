#pragma once

#include <cassert>
#include <vector>

#include <simdpp/simd.h>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

struct Row;

#if SIMDPP_USE_AVX2
constexpr int S = 64;
#elif SIMDPP_USE_SSE2
constexpr int S = 32;
#else
constexpr int S = 32;
#endif


namespace T
{

using uint = uint8_t;
using pack = simdpp::uint16<S / sizeof(uint), void>;
using mask = simdpp::mask_int16<S / sizeof(uint), void>;
using vec_pack = std::vector<pack, simdpp::aligned_allocator<pack, sizeof(pack)> >;
using matrix = std::vector<vec_pack>;
using vec_uint = std::vector<uint, simdpp::aligned_allocator<uint, sizeof(uint)> >;
using row = Row; // A row of vectors that can be run in parallel
using arr = std::array<row, 4>;
using arr_uint = std::array<uint, S / sizeof(uint)>;

} // namespace T


template<typename Tpack, typename Tuint>
inline Tpack
shift_one_right(Tpack pack, Tuint const left)
{
  #if SIMDPP_USE_AVX2
  std::vector<Tuint> vec(Tpack::length + 1, left);
  simdpp::store_u(&vec[1], pack);
  return simdpp::load_u(&vec[0]);
  #else
  return simdpp::align8<7, 8>(static_cast<Tpack>(simdpp::make_uint(left)), pack);
  #endif
}


template<typename Tpack>
inline Tpack
shift_one_right(Tpack pack)
{
  return simdpp::move8_r<1>(pack);
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


template <typename Trow>
void
print_score_vector_standard(Trow const & vX)
{
  using Tuint = typename Trow::Tuint;
  using Vec = std::vector<Tuint, simdpp::aligned_allocator<Tuint, sizeof(Tuint)> >;
  using Mat = std::vector<Vec>;

  long const t = vX.vectors.size();

  if (t == 0)
    return;

  Vec vec(vX.vectors[0].length, 0);
  Mat m(vX.vectors.size(), vec);

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


template <typename Tpack>
void
print_score_vectors(Tpack const & vH,
                    Tpack const & vH_up,
                    Tpack const & vE,
                    Tpack const & vF,
                    Tpack const & vF_up,
                    Tpack const & vW
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
