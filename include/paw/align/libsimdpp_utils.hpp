#pragma once

#include <array>
#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>

#include <simdpp/simd.h>


namespace paw
{

#if SIMDPP_USE_AVX512BW
constexpr int S = 64;
#elif SIMDPP_USE_AVX2
constexpr int S = 32;
#elif SIMDPP_USE_SSE2
constexpr int S = 16;
#else
constexpr int S = 16;
#endif

template <typename Tuint>
struct T : std::false_type
{};

template <>
struct T<uint8_t>
{
  using pack = simdpp::uint8<S / sizeof(uint8_t), void>;
  using mask = pack::mask_vector_type;
  using uint = pack::element_type;
  using vec_pack = std::vector<pack>;
  using vec_uint = std::vector<uint>;
  using arr_uint = std::array<uint, S / sizeof(uint)>;
  using arr_vec_pack = std::array<vec_pack, 5>;
};

template <>
struct T<uint16_t>
{
  using pack = simdpp::uint16<S / sizeof(uint16_t), void>;
  using mask = pack::mask_vector_type;
  using uint = pack::element_type;
  using vec_pack = std::vector<pack>;
  using vec_uint = std::vector<uint>;
  using arr_uint = std::array<uint, S / sizeof(uint)>;
  using arr_vec_pack = std::array<vec_pack, 5>;
};


namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tuint>
inline typename T<Tuint>::pack
shift_one_right(typename T<Tuint>::pack pack,
                typename T<Tuint>::uint const left,
                std::array<long, S / sizeof(typename T<Tuint>::uint)> const & reductions
                )
{
  std::array<typename T<Tuint>::uint, T<Tuint>::pack::length + 1> vec;
  vec[0] = left;
  simdpp::store_u(&vec[1], pack);
  Tuint const min_value = std::numeric_limits<Tuint>::min();

  for (long e = 1; e < static_cast<long>(T<Tuint>::pack::length); ++e)
  {
    long const val = static_cast<long>(vec[e]) + reductions[e - 1] - reductions[e];
    vec[e] = val >= min_value ? val : min_value;
  }

  return simdpp::load_u(&vec[0]);
}


inline long
magic_function(char const c)
{
  switch (c)
  {
  case 'n':
  case 'N':
    return 4;

  default:
    return 0x03 & ((c >> 2) ^ (c >> 1));
  }
}


template <typename Tuint>
inline void
init_vH_up(typename T<Tuint>::vec_pack & vH_up,
           Tuint const gap_open_val
           )
{
  typename T<Tuint>::vec_uint new_vH0(T<Tuint>::pack::length, 2 * gap_open_val + std::numeric_limits<Tuint>::min());
  assert(vH_up.size() > 0);
  new_vH0[0] = gap_open_val * 3 + std::numeric_limits<Tuint>::min();
  vH_up[0] = simdpp::load_u(&new_vH0[0]);
}


template <typename Tuint>
inline typename T<Tuint>::mask
max_greater(typename T<Tuint>::pack & v1, typename T<Tuint>::pack const & v2)
{
  typename T<Tuint>::mask is_greater = v2 > v1;
  v1 = simdpp::blend(v2, v1, is_greater);
  return is_greater;
}


template <typename Tuint>
void
print_pack(typename T<Tuint>::pack const & pack)
{
  using T = typename T<Tuint>::uint;

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


template <typename Tuint>
inline void
print_score_vector_standard(long m,
                            typename T<Tuint>::vec_pack const & vX,
                            long const top_left_score,
                            long const x_gain,
                            long const y_gain_total,
                            std::array<long, S / sizeof(Tuint)> const & reductions
                            )
{
  using Tvec_uint = typename T<Tuint>::vec_uint;

  long const t = vX.size();

  if (t == 0)
    return;

  Tvec_uint vec(vX[0].length, 0);
  std::vector<Tvec_uint> mat(t, vec);

  for (long v = 0; v < t; ++v)
  {
    simdpp::store_u(&mat[v][0], vX[v]);
  }

  for (long j = 0; j <= m; ++j)
  {
    long const v = j % t;
    long const e = j / t;
    assert(v < static_cast<long>(mat.size()));
    assert(e < static_cast<long>(mat[v].size()));

    if (v == 0 && j > 0)
      std::cout << "|" << std::setw(3);
    else
      std::cout << std::setw(4);

    std::cout << static_cast<long>(mat[v][e] - top_left_score - y_gain_total - x_gain * j + reductions[e]) << " ";
  }

  std::cout << "(" << t << " vectors, " << m << " elements)\n";
}


template <typename Tuint>
inline void
print_score_vectors(long m,
                    typename T<Tuint>::vec_pack const & vH,
                    typename T<Tuint>::vec_pack const & vH_up,
                    typename T<Tuint>::vec_pack const & vE,
                    typename T<Tuint>::vec_pack const & vF,
                    typename T<Tuint>::vec_pack const & vF_up,
                    typename T<Tuint>::vec_pack const & vW,
                    Tuint const top_left_score,
                    Tuint const x_gain,
                    Tuint const y_gain_total,
                    std::array<long, S / sizeof(Tuint)> const & reductions
                    )
{
  std::cout << "Standard H_up  : "; print_score_vector_standard<Tuint>(m,
                                                                       vH_up,
                                                                       top_left_score,
                                                                       x_gain,
                                                                       y_gain_total,
                                                                       reductions);
  std::cout << "Standard H     : "; print_score_vector_standard<Tuint>(m,
                                                                       vH,
                                                                       top_left_score,
                                                                       x_gain,
                                                                       y_gain_total,
                                                                       reductions);
  std::cout << "Standard E     : "; print_score_vector_standard<Tuint>(m,
                                                                       vE,
                                                                       top_left_score,
                                                                       x_gain,
                                                                       y_gain_total,
                                                                       reductions);
  std::cout << "Standard F_up  : "; print_score_vector_standard<Tuint>(m,
                                                                       vF_up,
                                                                       top_left_score,
                                                                       x_gain,
                                                                       y_gain_total,
                                                                       reductions);
  std::cout << "Standard F     : "; print_score_vector_standard<Tuint>(m,
                                                                       vF,
                                                                       top_left_score,
                                                                       x_gain,
                                                                       y_gain_total,
                                                                       reductions);
  std::cout << "Standard W     : "; print_score_vector_standard<Tuint>(m,
                                                                       vW,
                                                                       top_left_score,
                                                                       x_gain,
                                                                       y_gain_total,
                                                                       reductions);
  std::cout << "=====\n";
}


} //namespace SIMDPP_ARCH_NAMESPACE
} // anon namespace
