#pragma once

#include <array>
#include <cassert>
#include <vector>

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
  using uint = pack::uint_element_type;
  using vec_pack = std::vector<pack>;
  using vec_uint = std::vector<uint>;
  using arr_uint = std::array<uint, S / sizeof(uint)>;
  using arr_vec_pack = std::array<vec_pack, 4>;
};

template <>
struct T<uint16_t>
{
  using pack = simdpp::uint16<S / sizeof(uint16_t), void>;
  using mask = pack::mask_vector_type;
  using uint = pack::uint_element_type;
  using vec_pack = std::vector<pack>;
  using vec_uint = std::vector<uint>;
  using arr_uint = std::array<uint, S / sizeof(uint)>;
  using arr_vec_pack = std::array<vec_pack, 4>;
};


namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tuint>
inline typename T<Tuint>::pack
shift_one_right(typename T<Tuint>::pack pack,
                typename T<Tuint>::uint const left
                )
{
  std::array<typename T<Tuint>::uint, T<Tuint>::pack::length + 1> vec;
  //vec.fill(left);
  vec[0] = left;
  simdpp::store_u(&vec[1], pack);
  return simdpp::load_u(&vec[0]);
}


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

  for (long e = 1; e < static_cast<long>(T<Tuint>::pack::length); ++e)
  {
    long const val = static_cast<long>(vec[e]) + reductions[e - 1] - reductions[e];
    vec[e] = val > 0 ? val : 0;
  }

  return simdpp::load_u(&vec[0]);
}


template <typename Tuint>
inline typename T<Tuint>::pack
shift_one_right(typename T<Tuint>::pack pack)
{
  return shift_one_right(pack, 0);
}


inline long
magic_function(char const c)
{
  return 0x03 & ((c >> 2) ^ (c >> 1));
}


template <typename Tuint>
inline void
init_vH_up(typename T<Tuint>::vec_pack & vH_up,
           Tuint const gap_open_val
           )
{
  typename T<Tuint>::vec_uint new_vH0(T<Tuint>::pack::length, gap_open_val);
  assert(vH_up.size() > 0);
  new_vH0[0] = gap_open_val * 2;
  vH_up[0] = simdpp::load_u(&new_vH0[0]);
}


template <typename Tuint>
inline typename T<Tuint>::mask
max_greater(typename T<Tuint>::pack & v1, typename T<Tuint>::pack const & v2)
{
  // TODO: Only use simdpp::max when no backtracking.
  typename T<Tuint>::mask is_greater = v2 > v1;
  v1 = simdpp::blend(v2, v1, is_greater);
  return is_greater;
}


template <typename Tpack>
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
inline void
print_score_vector_standard(Trow const & vX)
{
  long const t = vX.vectors.size();

  if (t == 0)
    return;

  typename Trow::vec_uint vec(vX.vectors[0].length, 0);
  std::vector<typename Trow::vec_uint> m(vX.vectors.size(), vec);

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


template <typename Trow>
inline void
print_score_vectors(Trow const & vH,
                    Trow const & vH_up,
                    Trow const & vE,
                    Trow const & vF,
                    Trow const & vF_up,
                    Trow const & vW
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
