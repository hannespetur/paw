#pragma once

#include <cstdint>
#include <string> // std::string
#include <vector> // std::vector<T>
#include <iostream>

#include <simdpp/simd.h>

#include <paw/align/libsimdpp_utils.hpp>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

enum CigarOperation
{
  MATCH = 0,
  INSERTION,
  DELETION
};


struct Cigar
{
  std::size_t count;
  CigarOperation operation;
};


struct Backtrack
{
  std::size_t static const BACKTRACKS_PER_BYTE = 2;
  std::size_t static const BT_PER_CELL = sizeof(T::uint) * BACKTRACKS_PER_BYTE;

  /// \short How many bits are required for one column of backtracks
  std::size_t static const N_BT_BITS = 8 / BACKTRACKS_PER_BYTE;

  T::uint static const DEL_SHIFT = 0;
  T::uint static const INS_SHIFT = 1;
  T::uint static const DEL_E_SHIFT = 2;
  T::uint static const INS_E_SHIFT = 3;

  T::uint static const SUB_BT = 0; // substitution is represented with 0
  T::uint static const DEL_BT = 1 << DEL_SHIFT;
  T::uint static const INS_BT = 1 << INS_SHIFT;
  T::uint static const DEL_E_BT = 1 << DEL_E_SHIFT;
  T::uint static const INS_E_BT = 1 << INS_E_SHIFT;

  T::matrix matrix;
  T::arr_uint tmp_vec;

  Backtrack()
    : Backtrack(0, 0)
  {}


  Backtrack(std::size_t const n_row, std::size_t const n_vectors)
    : matrix(n_row, {(n_vectors + BT_PER_CELL - 1) / BT_PER_CELL, simdpp::make_zero()})
  {}

  void inline
  set_del(std::size_t const i /*row index*/,
          std::size_t const v /*vector index*/,
          T::mask const mask /*mask to set*/
          )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] |
        (simdpp::blend(static_cast<T::pack>(simdpp::make_uint(1 << (N_BT_BITS * (v % BT_PER_CELL) + DEL_SHIFT))),
                       static_cast<T::pack>(simdpp::make_zero()),
                       mask)
         );
  }


  void inline
  set_ins(std::size_t const i /*row index*/,
          std::size_t const v /*vector index*/,
          T::mask const mask /*mask to set*/
          )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] |
        (simdpp::blend(static_cast<T::pack>(simdpp::make_uint(1 << (N_BT_BITS * (v % BT_PER_CELL) + INS_SHIFT))),
                       static_cast<T::pack>(simdpp::make_zero()),
                       mask)
         );
  }


  void inline
  set_del_extend(std::size_t const i /*row index*/,
                 std::size_t const v /*vector index*/,
                 T::mask const mask /*mask to set*/
                 )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] |
        (simdpp::blend(static_cast<T::pack>(simdpp::make_uint(1 << (N_BT_BITS * (v % BT_PER_CELL) + DEL_E_SHIFT))),
                       static_cast<T::pack>(simdpp::make_zero()),
                       mask)
         );
  }


  void inline
  set_ins_extend(std::size_t const i /*row index*/,
                 std::size_t const v /*vector index*/,
                 T::mask const mask /*mask to set*/
                 )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] |
        (simdpp::blend(static_cast<T::pack>(simdpp::make_uint(1 << (N_BT_BITS * (v % BT_PER_CELL) + INS_E_SHIFT))),
                       static_cast<T::pack>(simdpp::make_zero()),
                       mask)
         );
  }


  bool inline
  is_del(std::size_t const i /*row index*/,
         std::size_t const v /*vector index*/,
         std::size_t const e /*element index*/
         )
  {
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return tmp_vec[e] & (DEL_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


  bool inline
  is_ins(std::size_t const i /*row index*/,
         std::size_t const v /*vector index*/,
         std::size_t const e /*element index*/
         )
  {
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return tmp_vec[e] & (INS_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


  bool inline
  is_del_extend(std::size_t const i /*row index*/,
                std::size_t const v /*vector index*/,
                std::size_t const e /*element index*/
                )
  {
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return tmp_vec[e] & (DEL_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


  bool inline
  is_ins_extend(std::size_t const i /*row index*/,
                std::size_t const v /*vector index*/,
                std::size_t const e /*element index*/
                )
  {
    simdpp::store_u(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    return tmp_vec[e] & (INS_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


};


template <typename Tint>
std::ostream &
operator<<(std::ostream & ss, std::vector<Cigar> const & cigar);


} // namespace SIMDPP_ARCH_NAMESPACE
} //namespace paw


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
    case MATCH:
      std::cout << "M"; break;

    case INSERTION:
      std::cout << "I"; break;

    case DELETION:
      std::cout << "D"; break;
    }
  }

  return ss;
}


std::size_t constexpr Backtrack::BACKTRACKS_PER_BYTE;
std::size_t constexpr Backtrack::BT_PER_CELL;
std::size_t constexpr Backtrack::N_BT_BITS;
T::uint constexpr Backtrack::DEL_SHIFT;
T::uint constexpr Backtrack::INS_SHIFT;
T::uint constexpr Backtrack::DEL_E_SHIFT;
T::uint constexpr Backtrack::INS_E_SHIFT;
T::uint constexpr Backtrack::SUB_BT;
T::uint constexpr Backtrack::DEL_BT;
T::uint constexpr Backtrack::INS_BT;
T::uint constexpr Backtrack::DEL_E_BT;
T::uint constexpr Backtrack::INS_E_BT;


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#endif // IMPLEMENT_PAW
