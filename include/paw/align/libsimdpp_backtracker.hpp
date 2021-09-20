#pragma once

#include <cstdint>
#include <string> // std::string
#include <vector> // std::vector<T>
#include <iostream>

#include <simdpp/simd.h>

#include <paw/align/cigar.hpp>
#include <paw/align/libsimdpp_utils.hpp>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tuint>
struct Backtrack
{
  using Tpack = typename T<Tuint>::pack;
  using Tmask = typename T<Tuint>::mask;
  using Tvec_pack = typename T<Tuint>::vec_pack;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  std::size_t static constexpr BACKTRACKS_PER_BYTE = 2;
  std::size_t static constexpr BT_PER_CELL = sizeof(Tuint) * BACKTRACKS_PER_BYTE;

  /// \short How many bits are required for one column of backtracks
  std::size_t static constexpr N_BT_BITS = 8 / BACKTRACKS_PER_BYTE;

  Tuint static constexpr DEL_SHIFT = 0;
  Tuint static constexpr INS_SHIFT = 1;
  Tuint static constexpr DEL_E_SHIFT = 2;
  Tuint static constexpr INS_E_SHIFT = 3;
  Tuint static constexpr SUB_BT = 0; // substitution is represented with 0
  Tuint static constexpr DEL_BT = 1 << DEL_SHIFT;
  Tuint static constexpr INS_BT = 1 << INS_SHIFT;
  Tuint static constexpr DEL_E_BT = 1 << DEL_E_SHIFT;
  Tuint static constexpr INS_E_BT = 1 << INS_E_SHIFT;

  long t{0};
  std::vector<Tvec_pack> matrix;
  mutable Tarr_uint tmp_vec;

  Backtrack()
    : Backtrack(0, 0)
  {}


  Backtrack(long const n_row, long const n_vectors)
    : t(n_vectors)
    , matrix(static_cast<std::size_t>(n_row),
             {(n_vectors + BT_PER_CELL - 1) / BT_PER_CELL, static_cast<Tpack>(simdpp::make_zero())})
  {
    assert(n_row >= 0);
    tmp_vec.fill(0);
  }


  inline Tpack const &
  get_pack(long const i /*row index*/, long const v /*vector index*/) const
  {
    assert(i < static_cast<long>(matrix.size()));
    return matrix[i][v / BT_PER_CELL];
  }


  void inline
  set_del(long const i /*row index*/,
          long const v /*vector index*/,
          Tmask const mask /*mask to set*/)
  {
    assert(i < static_cast<long>(matrix.size()));
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL]
                                 | (simdpp::blend(static_cast<Tpack>(simdpp::make_uint(DEL_BT <<
                                                                                       (N_BT_BITS *
                                                                                        (v % BT_PER_CELL)))),
                                                  static_cast<Tpack>(simdpp::make_zero()),
                                                  mask)
                                    );
  }


  void inline
  set_ins(long const i /*row index*/,
          long const v /*vector index*/,
          Tmask const mask /*mask to set*/)
  {
    assert(i < static_cast<long>(matrix.size()));
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] |
                                 (simdpp::blend(static_cast<Tpack>(simdpp::make_uint(INS_BT <<
                                                                                     (N_BT_BITS * (v % BT_PER_CELL)))),
                                                static_cast<Tpack>(simdpp::make_zero()),
                                                mask)
                                 );
  }


  void inline
  set_del_extend(long const i /*row index*/,
                 long const v /*vector index*/,
                 Tmask const mask /*mask to set*/)
  {
    assert(i < static_cast<long>(matrix.size()));
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] |
                                 (simdpp::blend(static_cast<Tpack>(simdpp::make_uint(DEL_E_BT <<
                                                                                     (N_BT_BITS * (v % BT_PER_CELL)))),
                                                static_cast<Tpack>(simdpp::make_zero()),
                                                mask)
                                 );
  }


  void inline
  set_ins_extend(long const i /*row index*/,
                 long const v /*vector index*/,
                 Tmask const mask /*mask to set*/)
  {
    assert(i < static_cast<long>(matrix.size()));
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] |
                                 (simdpp::blend(static_cast<Tpack>(simdpp::make_uint(INS_E_BT <<
                                                                                     (N_BT_BITS * (v % BT_PER_CELL)))),
                                                static_cast<Tpack>(simdpp::make_zero()),
                                                mask)
                                 );
  }


  /// \short Checks if an element is a deletion
  /// \param[in] i row index
  /// \param[in] v vector index
  /// \param[in] e element index
  /// \return True is the element is a deletion, otherwise false
  bool inline
  is_del(long i,
         long const v,
         long const e
         ) const
  {
    assert(i < static_cast<long>(matrix.size()));
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return (tmp_vec[e] >> (N_BT_BITS * (v % BT_PER_CELL))) & DEL_BT;
  }


  /// \short Checks if an element is an insertion
  /// \param[in] i row index
  /// \param[in] v vector index
  /// \param[in] e element index
  /// \return True is the element is an insertion, otherwise false
  bool inline
  is_ins(long const i,
         long const v,
         long const e
         ) const
  {
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return (tmp_vec[e] >> (N_BT_BITS * (v % BT_PER_CELL))) & INS_BT;
  }


  bool inline
  is_del_extend(long const i /*row index*/,
                long const v /*vector index*/,
                long const e /*element index*/
                ) const
  {
    assert(i < static_cast<long>(matrix.size()));
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return (tmp_vec[e] >> (N_BT_BITS * (v % BT_PER_CELL))) & DEL_E_BT;
  }


  bool inline
  is_ins_extend(long const i /*row index*/,
                long const v /*vector index*/,
                long const e /*element index*/
                ) const
  {
    assert(i < static_cast<long>(matrix.size()));
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return (tmp_vec[e] >> (N_BT_BITS * (v % BT_PER_CELL))) & INS_E_BT;
  }


};


template <typename Tint>
std::ostream &
operator<<(std::ostream & ss, std::vector<Cigar> const & cigar);


//inline void print_backtrack(Backtrack const & mB);


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#if defined(IMPLEMENT_PAW)

namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tint>
std::ostream &
operator<<(std::ostream & ss, std::vector<Cigar> const & cigar)
{
  for (auto const & c : cigar)
  {
    std::cout << c.count;

    switch (c.operation)
    {
    case CigarOperation::MATCH:
      std::cout << "M"; break;

    case CigarOperation::INSERTION:
      std::cout << "I"; break;

    case CigarOperation::DELETION:
      std::cout << "D"; break;

    default:
      assert(false);
    }
  }

  return ss;
}


/*
inline void
print_backtrack(Backtrack const & mB)
{
  std::cout << mB.BT_PER_CELL << "\n";
  std::cout << mB.matrix[0].size() << "\n";

  for (long i = 0; i < static_cast<long>(mB.matrix.size()); ++i)
  {
    long const t = mB.matrix[0].size();

    T::vec_uint vec(S / sizeof(T::uint), 0);
    std::vector<T::vec_uint> mat(t, vec);

    for (long v = 0; v < t; ++v)
    {
      T::pack const & pack = mB.get_pack(i, v);
      simdpp::store_u(&mat[v][0], pack);
    }

    for (long j = 0; j < t * T::pack::length; ++j)
    {
      std::size_t const v = j % t;
      std::size_t const e = j / t;
      T::uint const element = mat[v][e];

      //for (long k = 0; k < static_cast<long>(mB.BT_PER_CELL); ++k)
      //{
      //  std::cout << std::setw(1) << std::hex << static_cast<uint64_t>((element >> (4 * k)) & static_cast<T::uint>(0x000F));
      //}
      for (long k = mB.BT_PER_CELL - 1; k >= 0; --k)
      {
        std::cout << std::setw(1) << std::hex << static_cast<uint64_t>((element >> (4 * k)) & static_cast<T::uint>(0x000F));
      }

      std::cout << "|";
    }

    std::cout << std::endl;
  }
}
*/

template <typename Tuint>
std::size_t constexpr Backtrack<Tuint>::BACKTRACKS_PER_BYTE;
template <typename Tuint>
std::size_t constexpr Backtrack<Tuint>::BT_PER_CELL;
template <typename Tuint>
std::size_t constexpr Backtrack<Tuint>::N_BT_BITS;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::DEL_SHIFT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::INS_SHIFT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::DEL_E_SHIFT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::INS_E_SHIFT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::SUB_BT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::DEL_BT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::INS_BT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::DEL_E_BT;
template <typename Tuint>
Tuint constexpr Backtrack<Tuint>::INS_E_BT;


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#endif // IMPLEMENT_PAW
