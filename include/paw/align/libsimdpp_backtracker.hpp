#pragma once

#include <cstdint>
#include <string> // std::string
#include <vector> // std::vector<T>

#include "event.hpp"

//#include <boost/simd.hpp>
//#include <boost/simd/memory/allocator.hpp>

#include <simdpp/simd.h>


namespace paw
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


template <typename Tuint>
struct Backtrack
{
  using Tint = Tuint;
  using Tpack = simdpp::int16<16 / sizeof(Tint), void>;
  using Tvec = std::vector<Tpack, simdpp::aligned_allocator<Tpack, sizeof(Tpack)> >;
  using Tmatrix = std::vector<Tvec>;

  std::size_t static const BACKTRACKS_PER_BYTE = 2;
  std::size_t static const BT_PER_CELL = sizeof(Tuint) * BACKTRACKS_PER_BYTE;
  std::size_t static const N_BT_BITS = 8 / BACKTRACKS_PER_BYTE;

  Tuint static const DEL_SHIFT = 0;
  Tuint static const INS_SHIFT = 1;
  Tuint static const DEL_E_SHIFT = 2;
  Tuint static const INS_E_SHIFT = 3;

  Tuint static const SUB_BT = 0; // substitution is represented with 0
  Tuint static const DEL_BT = 1 << DEL_SHIFT;
  Tuint static const INS_BT = 1 << INS_SHIFT;
  Tuint static const DEL_E_BT = 1 << DEL_E_SHIFT;
  Tuint static const INS_E_BT = 1 << INS_E_SHIFT;

  Tmatrix matrix;

  Backtrack()
    : Backtrack(0, 0)
  {}


  Backtrack(std::size_t const n_row, std::size_t const n_vectors)
    : matrix(n_row, {(n_vectors + BT_PER_CELL - 1) / BT_PER_CELL, simdpp::make_zero()})
  {}

  void inline
  set_del(std::size_t const i /*row index*/,
          std::size_t const v /*vector index*/,
          Tpack const pack /*pack to set*/
          )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] | (pack << (N_BT_BITS * (v % BT_PER_CELL) + DEL_SHIFT));
  }


  void inline
  set_ins(std::size_t const i /*row index*/,
          std::size_t const v /*vector index*/,
          Tpack const pack /*pack to set*/
          )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] | (pack << (N_BT_BITS * (v % BT_PER_CELL) + INS_SHIFT));
  }


  void inline
  set_del_extend(std::size_t const i /*row index*/,
                 std::size_t const v /*vector index*/,
                 Tpack const pack /*pack to set*/
                 )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] | (pack << (N_BT_BITS * (v % BT_PER_CELL) + DEL_E_SHIFT));
  }


  void inline
  set_ins_extend(std::size_t const i /*row index*/,
                 std::size_t const v /*vector index*/,
                 Tpack const pack /*pack to set*/
                 )
  {
    matrix[i][v / BT_PER_CELL] = matrix[i][v / BT_PER_CELL] | (pack << (N_BT_BITS * (v % BT_PER_CELL) + INS_E_SHIFT));
  }


  bool inline
  is_del(std::size_t const i /*row index*/,
         std::size_t const v /*vector index*/,
         std::size_t const e /*element index*/
         ) const
  {
    //Tpack tmp = matrix[i][v / BT_PER_CELL] & (DEL_BT << (N_BT_BITS * (v % BT_PER_CELL)));
    std::vector<Tuint, simdpp::aligned_allocator<Tuint, sizeof(Tuint)> > tmp_vec(16 / sizeof(Tint));
    simdpp::store(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    //
    //std::cerr << "IS DEL? e = " << e << "\n";
    //
    //for (auto const item : tmp_vec)
    //  std::cerr << item << ",  ";
    //
    //std::cerr << std::endl;
    //
    //assert(e < tmp_vec.size());
    return tmp_vec[e] & (DEL_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


  bool inline
  is_ins(std::size_t const i /*row index*/,
         std::size_t const v /*vector index*/,
         std::size_t const e /*element index*/
         ) const
  {
    //Tpack tmp = matrix[i][v / BT_PER_CELL] & (INS_BT << (N_BT_BITS * (v % BT_PER_CELL)));
    std::vector<Tuint, simdpp::aligned_allocator<Tuint, sizeof(Tuint)> > tmp_vec(16 / sizeof(Tint));
    simdpp::store(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);

    //std::cerr << "IS INS? e = " << e << "\n";
    //
    //for (auto const item : tmp_vec)
    //  std::cerr << item << ",  ";
    //
    //std::cerr << std::endl;

    assert(e < tmp_vec.size());
    return tmp_vec[e] & (INS_BT << (N_BT_BITS * (v % BT_PER_CELL)));
    //return matrix[i][v / BT_PER_CELL][e] & (INS_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


  bool inline
  is_del_extend(std::size_t const i /*row index*/,
                std::size_t const v /*vector index*/,
                std::size_t const e /*element index*/
                ) const
  {
    //Tpack tmp = matrix[i][v / BT_PER_CELL] & (DEL_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
    std::vector<Tuint, simdpp::aligned_allocator<Tuint, sizeof(Tuint)> > tmp_vec(16 / sizeof(Tint));
    simdpp::store(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);

    //std::cerr << "e = " << e << "\n";
    //
    //for (auto const item : tmp_vec)
    //  std::cerr << item << ",  ";
    //
    //std::cerr << std::endl;
    //
    //assert(e < tmp_vec.size());

    return tmp_vec[e] & (DEL_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
    //return matrix[i][v / BT_PER_CELL][e] & (DEL_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


  bool inline
  is_ins_extend(std::size_t const i /*row index*/,
                std::size_t const v /*vector index*/,
                std::size_t const e /*element index*/
                ) const
  {
    //Tpack tmp = matrix[i][v / BT_PER_CELL] & (INS_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
    std::vector<Tuint, simdpp::aligned_allocator<Tuint, sizeof(Tuint)> > tmp_vec(16 / sizeof(Tint));
    simdpp::store(&tmp_vec[0], matrix[i][v / BT_PER_CELL]);
    //assert(e < tmp_vec.size());
    return tmp_vec[e] & (INS_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
    //return matrix[i][v / BT_PER_CELL][e] & (INS_E_BT << (N_BT_BITS * (v % BT_PER_CELL)));
  }


};


template <typename Tint>
std::ostream &
operator<<(std::ostream & ss, std::vector<Cigar> const & cigar);


} //namespace paw


#if defined(IMPLEMENT_PAW)


namespace paw
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


//std::size_t static constexpr BACKTRACKS_PER_BYTE = 2;
//std::size_t static constexpr BT_PER_CELL = sizeof(Tuint) * BACKTRACKS_PER_BYTE;
//std::size_t static constexpr N_BT_BITS = 8 / BACKTRACKS_PER_BYTE;

template <typename Tuint>
std::size_t constexpr Backtracker<Tuint>::BACKTRACKS_PER_BYTE;

template <typename Tuint>
std::size_t constexpr Backtracker<Tuint>::BT_PER_CELL;

template <typename Tuint>
std::size_t constexpr Backtracker<Tuint>::N_BT_BITS;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::DEL_SHIFT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::INS_SHIFT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::DEL_E_SHIFT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::INS_E_SHIFT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::SUB_BT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::DEL_BT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::INS_BT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::DEL_E_BT;

template <typename Tuint>
Tuint constexpr Backtracker<Tuint>::INS_E_BT;


} //namespace paw


#endif // IMPLEMENT_PAW
