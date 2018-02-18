#pragma once

#include <array> // std::array
#include <cstdint> // uint8_t
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <boost/simd/pack.hpp>
#include <boost/simd/memory/allocator.hpp>
#include <boost/simd/function/all.hpp>
#include <boost/simd/function/if_zero_else.hpp>
#include <boost/simd/function/aligned_load.hpp>
#include <boost/simd/function/aligned_store.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/store.hpp>
#include <boost/simd/function/maximum.hpp>

#include <paw/align/aligner_options.hpp>
#include <paw/align/backtracker.hpp>
#include <paw/align/row.hpp>


namespace paw
{


template <typename Tuint, typename Tit>
class Aligner
{
public:
  using Trow = Row<Tuint>;   // A row of vectors that can be run in parallel
  using Tarr = std::array<Trow, 4>;
  using Tpack = typename Row<Tuint>::Tpack;
  using Tlogical_pack = typename Row<Tuint>::Tlogical_pack;

  // p is the length (or cardinality) of the SIMD vectors
  std::size_t static constexpr p = boost::simd::cardinal_of<Tpack>();


private:
  AlignerOptions<Tuint> const opt;
  std::vector<Event> free_snp_edits;
  std::vector<Event> free_del_edits;
  std::vector<Event> free_ins_edits;
  Tit d_begin;
  Tit d_end;
  Tit q_begin;
  Tit q_end;

  std::size_t const m;
  std::size_t alignment_end = m;
  std::size_t n = 0;
  std::size_t const t = 0;
  Tuint x_gain;
  Tuint y_gain;
  Tuint const gap_open_val_x;
  Tpack const gap_open_pack_x;
  Tuint const gap_open_val_y;
  Tpack const gap_open_pack_y;
  Tuint const gap_open_val;
  Tpack const gap_open_pack;
  Tuint const max_score_val;
  Tpack const max_score_pack;
  int64_t total_reductions = 0;
  Trow vH_up; // Previous H row
  Trow vH;    // Current H row
  Trow vE;    // Current E row
  Trow vF_up; // Previous F row
  Trow vF;    // Current F row
  Tarr W_profile;
  Backtracker<Tuint> mB; //(n /*n_row*/, t /*n_vectors in each score row*/);
  Tuint top_left_score = 0;

  Tpack max_greater(Tpack & v1, Tpack const & v2);
  // Tpack max_greater_or_equal(Tpack & v1, Tpack const & v2);
  void calculate_DNA_W_profile();
  void calculate_scores();

  // Same as 'calculate_scores', but assumes default options (less runtime checks)
  // void calculate_scores_default();
  void check_gap_extend_deletions();
  void check_gap_extend_deletions_with_backtracking(std::size_t const i);
  void init_score_vectors();


public:
  Aligner(Tit d_begin,
          Tit d_end,
          AlignerOptions<Tuint> const & _opt = AlignerOptions<Tuint>(true /*default options*/),
          std::set<Event> const & free_edits = std::set<Event>()
          );

  int64_t align(Tit q_begin, Tit q_end);   // Align query
  std::pair<std::string, std::string> get_aligned_strings();
  std::set<Event> get_edit_script(std::pair<std::string, std::string> const & s);
  std::vector<Cigar> get_cigar(std::pair<std::string, std::string> const & s);

};


} // namespace paw


/******************
 * IMPLEMENTATION *
 ******************/

#if defined(IMPLEMENT_PAW) || defined(__JETBRAINS_IDE__)

#include <cassert>
#include <iomanip>

namespace
{

template <typename Tuint>
inline paw::Row<Tuint> const &
get_W(const char b, std::array<paw::Row<Tuint>, 4> const & W_profile)
{
  switch (b)
  {
  case 'A':
    return W_profile[0];

  case 'C':
    return W_profile[1];

  case 'G':
    return W_profile[2];

  case 'T':
    return W_profile[3];

  default:
    std::cerr << "Warning: " << b << " is not A, C, G, or T.\n";
    assert(false);
    return W_profile[3];
  }
}


template <typename Tpack>
inline void
print_score_vector_standard(Tpack const & vX)
{
  std::size_t const t = vX.vectors.size();

  for (std::size_t j = 0; j < vX.n_elements; ++j)
  {
    std::size_t const v = j % t;
    std::size_t const e = j / t;
    std::cout << std::setw(4) << static_cast<uint32_t>(vX.vectors.at(v)[e]) << " ";
  }

  std::cout << "\n";
}


template <typename Tpack>
inline void
print_score_vector_vectorized(Tpack const & v)
{
  for (auto const & vec : v.vectors)
    std::cout << vec << " ";

  std::cout << "\n";
}


/*
template <typename Tpack>
inline void
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
  std::cout << "\n";
  std::cout << "Vectorized H_up: "; print_score_vector_vectorized(vH_up);
  std::cout << "Vectorized H   : "; print_score_vector_vectorized(vH);
  std::cout << "Vectorized E   : "; print_score_vector_vectorized(vE);
  std::cout << "=====\n";
}
*/


/*
constexpr int inline
shift_elements_left(int i, int c)
{
  return (c - 1 == i) ? -1 : (i + 1);
}
*/


constexpr int inline
shift_elements_right(int i, int /*c*/)
{
  return i - 1;
}


} // anon namespace


namespace paw
{

template <typename Tuint, typename Tit>
inline
Aligner<Tuint, Tit>::Aligner(Tit _d_begin /*database begin*/,
                             Tit _d_end /*database end*/,
                             AlignerOptions<Tuint> const & _opt,
                             std::set<Event> const & free_edits
                             )
  : opt(_opt)
  , d_begin(_d_begin)
  , d_end(_d_end)
  , m(std::distance(d_begin, d_end))
  , n(0)
  , t((m + p) / p)
  , x_gain(opt.gap_extend)
  , y_gain(std::max(opt.gap_extend, static_cast<Tuint>(opt.mismatch - x_gain)))
  , gap_open_val_x(opt.gap_open - x_gain)
  , gap_open_pack_x(gap_open_val_x)
  , gap_open_val_y(opt.gap_open - y_gain)
  , gap_open_pack_y(gap_open_val_y)
  , gap_open_val(std::max(gap_open_val_x, gap_open_val_y))
  , gap_open_pack(gap_open_val)
  , max_score_val(std::numeric_limits<Tuint>::max() - gap_open_val)
  , max_score_pack(max_score_val)
  , vH_up(m + 1)
  , vH(m + 1)
  , vE(m + 1)
  , vF_up(m + 1)
  , vF(m + 1)
  , W_profile
  {
    {
      Trow(m + 1, 0u),
      Trow(m + 1, 0u),
      Trow(m + 1, 0u),
      Trow(m + 1, 0u)
    }


  }
  , mB()
{
  for (auto const & e : free_edits)
  {
    if (e.ref.size() == 0)
    {
      free_ins_edits.push_back(e);
    }
    else if (e.alt.size() == 0)
    {
      free_del_edits.push_back(e);
    }
    else   //if (e.ref.size() == 1 && e.alt.size() == 1)
    {
      free_snp_edits.push_back(e);
    }
  }

  auto get_edit_str = [](Event const & e) -> std::string
                      {
                        std::stringstream ss;
                        ss << e.pos << " "
                           << (e.ref.size() > 0 ? e.ref : "-") << " "
                           << (e.alt.size() > 0 ? e.alt : "-");
                        return ss.str();
                      };


  if (free_snp_edits.size() > 0)
  {
    std::cout << "Free SNP edits are:\n";

    for (auto const & e : free_snp_edits)
      std::cout << get_edit_str(e) << "\n";
  }

  if (free_ins_edits.size() > 0)
  {
    std::cout << "Free INS edits are:\n";

    for (auto const & e : free_ins_edits)
      std::cout << get_edit_str(e) << "\n";
  }

  if (free_del_edits.size() > 0)
  {
    std::cout << "Free DEL edits are:\n";

    for (auto const & e : free_del_edits)
      std::cout << get_edit_str(e) << "\n";
  }

  calculate_DNA_W_profile();
}


/// Calculate the maximum with backtrack
template <typename Tuint, typename Tit>
typename Aligner<Tuint, Tit>::Tpack inline
Aligner<Tuint, Tit>::max_greater(Tpack & v1, Tpack const & v2)
{
  Tlogical_pack is_greater = v2 > v1;
  v1 = boost::simd::if_else(is_greater, v2, v1);
  return boost::simd::if_one_else_zero(is_greater);
}


/*
template <typename Tuint, typename Tit>
typename Aligner<Tuint, Tit>::Tpack inline
Aligner<Tuint, Tit>::max_greater_or_equal(Tpack & v1, Tpack const & v2)
{
  Tlogical_pack is_greater = v2 >= v1;
  v1 = boost::simd::if_else(is_greater, v2, v1);
  return boost::simd::if_one_else_zero(is_greater);
}
*/


/// Calculate W_profile
template <typename Tuint, typename Tit>
void inline
Aligner<Tuint, Tit>::calculate_DNA_W_profile()
{
  std::array<char, 4> constexpr DNA_BASES = {{'A', 'C', 'G', 'T'}};
  Tuint const match = opt.match;
  Tuint const mismatch = -opt.mismatch;
  std::size_t const t = W_profile[0].vectors.size();

  for (std::size_t i = 0; i < DNA_BASES.size(); ++i)
  {
    char const a = DNA_BASES[i];
    auto & W = W_profile[i];

    for (std::size_t v = 0; v < t; ++v)
    {
      for (std::size_t e = 0, j = v; j < m; ++e, j += t)
        W.vectors[v][e] = x_gain + y_gain + (a == *std::next(d_begin, j) ? match : mismatch);
    }

    // Update the W_profile if there are any free mismatches
    for (auto const & snp_e : free_snp_edits)
    {
      assert(snp_e.ref.size() == 1);
      assert(snp_e.alt.size() == 1);
      assert(snp_e.pos < m);
      assert(snp_e.ref[0] == *std::next(d_begin, snp_e.pos));

      if (a != snp_e.alt[0])
        continue;

      std::size_t const v = snp_e.pos % t;
      std::size_t const e = snp_e.pos / t;

      W.vectors[v][e] = x_gain + y_gain + match;
    }
  }
}


template <typename Tuint, typename Tit>
void inline
Aligner<Tuint, Tit>::init_score_vectors()
{
  if (opt.top_row_gap_open_free)
  {
    for (auto & v : vH_up.vectors)
      v = gap_open_val * 2;
  }
  else
  {
    for (auto & v : vH_up.vectors)
      v = gap_open_val;
  }

  vH.vectors = vH_up.vectors;

  for (auto & v : vF_up.vectors)
    v = 0;

  if (!opt.top_row_gap_open_free)
    vH_up.vectors[0][0] = gap_open_val * 2; // Set the very first value

  if (opt.top_row_free)
  {
    for (std::size_t v = 0; v < t; ++v)
    {
      for (std::size_t e = 0, j = v; j < m; ++e, j += t)
      {
        Tuint const gain = j * x_gain;
        vH_up.vectors[v][e] = 2 * gap_open_val + gain;
        vH.vectors[v][e] = 2 * gap_open_val + gain;
        vF_up.vectors[v][e] = gap_open_val + gain;
      }
    }
  }

  top_left_score = vH_up.vectors[0][0];
}


template <typename Tuint, typename Tit>
int64_t inline
Aligner<Tuint, Tit>::align(Tit _q_begin, Tit _q_end)
{
  q_begin = _q_begin;
  q_end = _q_end;
  n = std::distance(q_begin, q_end);
  init_score_vectors();
  total_reductions = 0;

  if (opt.backtracking)
  {
    if (n > mB.matrix.size())
    {
      mB = Backtracker<Tuint>(n, t);
    }
    else
    {
      // Make sure the first n rows have only zeros
      for (std::size_t i = 0; i < n; ++i)
      {
        for (auto & v : mB.matrix[i])
          v = 0x0;
      }
    }
  }

//  if (opt.default_options)
//    calculate_scores_default();
//  else
  calculate_scores();

  // Get final score
  if (opt.bottom_row_free)
  {
    std::size_t max_j = 0;
    int64_t max_s = std::numeric_limits<int64_t>::min();

    // Find the highest score on the bottom row
    for (std::size_t j = 0; j <= m; ++j)
    {
      int64_t s = static_cast<int64_t>(vH_up.vectors[j % t][j / t]) - j * x_gain;

      if (s >= max_s)
      {
        max_s = s;
        max_j = j;
      }
    }

    alignment_end = max_j;

    return static_cast<int64_t>(vH_up.vectors[max_j % t][max_j / t])
           + total_reductions
           - top_left_score
           - n * y_gain
           - max_j * x_gain;

  }
  else if (opt.bottom_row_gap_open_free)
  {
    std::size_t max_j = 0;
    int64_t max_s = std::numeric_limits<int64_t>::min();

    // Find the highest score on the bottom row
    for (std::size_t j = 0; j <= m; ++j)
    {
      int64_t s = static_cast<int64_t>(vH_up.vectors[j % t][j / t]);

      if (s >= max_s)
      {
        max_s = s;
        max_j = j;
      }
    }

    alignment_end = max_j;

    return static_cast<int64_t>(vH_up.vectors[max_j % t][max_j / t])
           + total_reductions
           - top_left_score
           - n * y_gain;
  }
  else
  {
    alignment_end = m;

    return static_cast<int64_t>(vH_up.vectors[m % t][m / t])
           + total_reductions
           - top_left_score
           - n * y_gain
           - m * x_gain;
  }
}


template <typename Tuint, typename Tit>
void inline
Aligner<Tuint, Tit>::check_gap_extend_deletions()
{
  for (std::size_t c = 0; c < 2; ++c)
  {
    Tuint max_val = 0;

    // Check for deletions in vector 0
    for (std::size_t e = 1; e < p; ++e)
    {
      max_val = std::max(max_val, static_cast<Tuint>(vE.vectors[t - 1][e - 1]));

      if (max_val > static_cast<Tuint>(vE.vectors[0][e]))
        vE.vectors[0][e] = max_val;
    }

    // Check for deletions in vectors 1,...,t-1
    for (std::size_t v = 1; v < t; ++v)
      vE.vectors[v] = boost::simd::max(vE.vectors[v - 1], vE.vectors[v]);
  }
}


template <typename Tuint, typename Tit>
void inline
Aligner<Tuint, Tit>::check_gap_extend_deletions_with_backtracking(std::size_t const i)
{
  Tpack vE0 {
    0
  };

  for (std::size_t c = 0; c < 2; ++c)
  {
    Tuint max_val = 0;

    // Check for deletions in vector 0
    for (std::size_t e = 1; e < p; ++e)
    {
      max_val = std::max(max_val, static_cast<Tuint>(vE.vectors[t - 1][e - 1]));

      if (max_val > static_cast<Tuint>(vE.vectors[0][e]))
        vE0[e] = max_val;
    }

    mB.set_del_extend(i, 0, max_greater(vE.vectors[0], vE0));

    // Check for deletions in vectors 1,...,t-1
    for (std::size_t v = 1; v < t; ++v)
    {
      mB.set_del_extend(i, v, max_greater(vE.vectors[v], vE.vectors[v - 1]));
    }
  }
}


template <typename Tuint, typename Tit>
void inline
Aligner<Tuint, Tit>::calculate_scores()
{
  bool const backtracking = opt.backtracking;
  bool const left_column_gap_open_free = opt.left_column_gap_open_free;
  bool const right_column_gap_open_free = opt.right_column_gap_open_free;
  std::size_t const right_v = m % t; // Vector that contains the rightmost element
  std::size_t const right_e = m / t; // The right-most element (in vector 'right_v')

  Tuint const max_gain_per_row = std::max(static_cast<Tuint>(0),
                                          static_cast<Tuint>(opt.match + x_gain + y_gain)
                                          );
  Tuint current_max_score = vH_up.vectors[0][0]; // Keep track of current max score

  /// Start of outer loop
  for (std::size_t i = 0; i < n; ++i)
  {
    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Trow const & vW = get_W(*std::next(q_begin, i), W_profile);

    // Clear vE vectors
    for (auto & vE_vec : vE.vectors)
      vE_vec = 0;

    /// Calculate vector 0
    vH.vectors[0] = boost::simd::shuffle<boost::simd::pattern<shift_elements_right> >(
      vH_up.vectors[t - 1] + vW.vectors[t - 1]);

    // This is the only way to reach the left-most coordinate
    vH.vectors[0][0] = std::max(static_cast<Tuint>(vF_up.vectors[0][0]),
                                static_cast<Tuint>(vH_up.vectors[0][0] - gap_open_val)
                                );

    // Check if any insertion have highest values
    if (backtracking)
    {
      mB.matrix[i][0][0] = mB.INS_BT;
      vF.vectors[0] = vH_up.vectors[0] - gap_open_pack;

      if (left_column_gap_open_free)
        vF.vectors[0][0] = vH_up.vectors[0][0];

      if (right_column_gap_open_free && right_v == 0)
        vF.vectors[0][right_e] = vH_up.vectors[0][right_e];

      mB.set_ins_extend(i, 0, max_greater(vF.vectors[0], vF_up.vectors[0]));
      mB.set_ins(i, 0, max_greater(vH.vectors[0], vF.vectors[0]));
    }
    else
    {
      vF.vectors[0] = vH_up.vectors[0] - gap_open_pack;

      if (left_column_gap_open_free)
        vF.vectors[0][0] = vH_up.vectors[0][0];

      if (right_column_gap_open_free && right_v == 0)
        vF.vectors[0][right_e] = vH_up.vectors[0][right_e];

      vF.vectors[0] = boost::simd::max(vF.vectors[0], vF_up.vectors[0]);
      vH.vectors[0] = boost::simd::max(vH.vectors[0], vF.vectors[0]);
    }
    /// Done calculating vector 0

    /// Calculate vectors 1,...,v-1
    for (std::size_t v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH.vectors[v] = vH_up.vectors[v - 1] + vW.vectors[v - 1];
      //vF.vectors[v] = boost::simd::max(vH_up.vectors[v] - gap_open_pack, vF_up.vectors[v]);

      if (backtracking)
      {
        vF.vectors[v] = vH_up.vectors[v] - gap_open_pack;

        // Slight performance gain if this if statement is moved out of the loop
        if (right_column_gap_open_free && v == right_v)
          vF.vectors[right_v][right_e] = vH_up.vectors[right_v][right_e];

        mB.set_ins_extend(i, v, max_greater(vF.vectors[v], vF_up.vectors[v]));
        mB.set_ins(i, v, max_greater(vH.vectors[v], vF.vectors[v]));
      }
      else
      {
        vF.vectors[v] = boost::simd::max(vH_up.vectors[v] - gap_open_pack, vF_up.vectors[v]);

        // Slight performance gain if this if statement is moved out of the loop
        if (right_column_gap_open_free && v == right_v)
          vF.vectors[right_v][right_e] = vH_up.vectors[right_v][right_e];

        vH.vectors[v] = boost::simd::max(vH.vectors[v], vF.vectors[v]);
      }

      // Deletions pass 1: Gap opens
      vE.vectors[v] = vH.vectors[v - 1] - gap_open_pack;
    } /// Done calculating vectors 1,...,v-1

    // Calculate vE.vector[0]
    vE.vectors[0] =
      boost::simd::shuffle<boost::simd::pattern<shift_elements_right> >(
        vH.vectors[t - 1] - gap_open_pack_x);

    //for (std::size_t e = 1; e < p; ++e)
    //  vE.vectors[0][e] = vH.vectors[t - 1][e - 1] - gap_open_val_x;

    // Deletions pass 2: Gap extends
    if (backtracking)
    {
      check_gap_extend_deletions_with_backtracking(i);

      // Check if any vE has higher scores than vH
      for (std::size_t v = 0; v < t; ++v)
        mB.set_del(i, v, max_greater(vH.vectors[v], vE.vectors[v]));
    }
    else
    {
      check_gap_extend_deletions();

      // Check if any vE has higher scores than vH
      for (std::size_t v = 0; v < t; ++v)
        vH.vectors[v] = boost::simd::max(vH.vectors[v], vE.vectors[v]);
    }

    //std::cout << "q[" << i << "] = " << q[i] << "\n";
    //print_score_vectors(vH, vH_up, vE, vF, vF_up, vW); // Useful when debugging
    //std::cout << "max score = " << static_cast<uint16_t>(current_max_score) << "\n";

    if (current_max_score + max_gain_per_row < max_score_val)
    {
      // Assume the worst-case and only find the max score when it would be possible to overflow
      current_max_score += max_gain_per_row;
    }
    else
    {
      // Find the correct max score (accurate but computationally expensive)
      current_max_score = 0;

      for (auto const & vec : vH.vectors)
      {
        current_max_score = std::max(current_max_score, boost::simd::maximum(vec));

        if (current_max_score >= max_score_val)
        {
          current_max_score -= gap_open_val;
          total_reductions += gap_open_val;

          for (auto & vH_vec : vH.vectors)
            vH_vec = boost::simd::max(vH_vec - gap_open_pack, gap_open_pack);

          for (auto & vF_vec : vF.vectors)
            vF_vec = boost::simd::max(vF_vec - gap_open_pack, gap_open_pack);

          break;
        }
      }
    }

    std::swap(vF.vectors, vF_up.vectors);
    std::swap(vH.vectors, vH_up.vectors);
  } /// End of outer loop
}

/*
template <typename Tuint, typename Tit>
void inline
Aligner<Tuint, Tit>::calculate_scores_default()
{
  std::size_t const right_v = m % t; // Vector that contains the rightmost element
  std::size_t const right_e = m / t; // The right-most element (in vector 'right_v')
  Tuint const max_gain_per_row = std::max(static_cast<Tuint>(0),
                                          static_cast<Tuint>(opt.match + x_gain + y_gain)
                                          );
  Tuint current_max_score = vH_up.vectors[0][0]; // Keep track of current max score

  /// Start of outer loop
  for (std::size_t i = 0; i < n; ++i)
  {
    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Trow const & vW = get_W(*std::next(q_begin, i), W_profile);

    // Clear vE vectors
    for (auto & vE_vec : vE.vectors)
      vE_vec = 0;

    /// Calculate vector 0
    // This is the only way to reach the left-most coordinate
    vH.vectors[0][0] = std::max(static_cast<Tuint>(vF_up.vectors[0][0]),
                                static_cast<Tuint>(vH_up.vectors[0][0] - gap_open_val)
                                );

    // Scalar check for substitutions and gap open deletions
    for (std::size_t e = 1; e < p; ++e)
      vH.vectors[0][e] = vH_up.vectors[t - 1][e - 1] + vW.vectors[t - 1][e - 1];

    // Check if any insertion have highest values
    mB.matrix[i][0][0] = mB.INS_BT;
    vF.vectors[0] = vH_up.vectors[0] - gap_open_pack;
    vF.vectors[0][0] = vH_up.vectors[0][0]; // left-most value is gap open free

    if (0 == right_v)
      vF.vectors[0][right_e] = vH_up.vectors[0][right_e];

    mB.set_ins_extend(i, 0, max_greater(vF.vectors[0], vF_up.vectors[0]));
    mB.set_ins(i, 0, max_greater(vH.vectors[0], vF.vectors[0]));
    /// Done calculating vector 0

    /// Calculate vectors 1,...,v-1
    for (std::size_t v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH.vectors[v] = vH_up.vectors[v - 1] + vW.vectors[v - 1];
      vF.vectors[v] = vH_up.vectors[v] - gap_open_pack;

      // Slight performance gain if this if statement is moved out of the loop
      if (v == right_v)
        vF.vectors[right_v][right_e] = vH_up.vectors[right_v][right_e];

      mB.set_ins_extend(i, v, max_greater(vF.vectors[v], vF_up.vectors[v]));
      mB.set_ins(i, v, max_greater(vH.vectors[v], vF.vectors[v]));

      // Deletions pass 1: Gap opens
      vE.vectors[v] = vH.vectors[v - 1] - gap_open_pack;
    } /// Done calculating vectors 1,...,v-1

    // Calculate vE.vector[0]
    for (std::size_t e = 1; e < p; ++e)
      vE.vectors[0][e] = vH.vectors[t - 1][e - 1] - gap_open_val_x;

    // Deletions pass 2: Gap extends
    check_gap_extend_deletions_with_backtracking(i);

    // Check if any vE has higher scores than vH
    for (std::size_t v = 0; v < t; ++v)
      mB.set_del(i, v, max_greater(vH.vectors[v], vE.vectors[v]));

    //std::cout << "q[" << i << "] = " << q[i] << "\n";
    //print_score_vectors(vH, vH_up, vE, vF, vF_up, vW); // Useful when debugging

    if (current_max_score + max_gain_per_row < max_score_val)
    {
      // Assume the worst-case and only find the max score when it would be possible to overflow
      current_max_score += max_gain_per_row;
    }
    else
    {
      // Find the correct max score (accurate but computationally expensive)
      current_max_score = 0;

      for (auto const & vec : vH.vectors)
      {
        current_max_score = std::max(current_max_score, boost::simd::maximum(vec));

        if (current_max_score >= max_score_val)
        {
          current_max_score -= gap_open_val;
          total_reductions += gap_open_val;

          for (auto & vH_vec : vH.vectors)
            vH_vec = boost::simd::max(vH_vec - gap_open_pack, gap_open_pack);

          for (auto & vF_vec : vF.vectors)
            vF_vec = boost::simd::max(vF_vec - gap_open_pack, gap_open_pack);

          break;
        }
      }
    }

    std::swap(vF.vectors, vF_up.vectors);
    std::swap(vH.vectors, vH_up.vectors);
  } /// End of outer loop
}
*/


template <typename Tuint, typename Tit>
std::pair<std::string, std::string> inline
Aligner<Tuint, Tit>::get_aligned_strings()
{
  std::size_t i = n;
  std::size_t j = alignment_end; //d.size();
  std::pair<std::string, std::string> s;

  // TODO: Fix this mess
  std::string d(d_begin, d_end);
  std::string q(q_begin, q_end);

  if (alignment_end < d.size())
  {
    s.first.append(std::string(d.rbegin(), d.rbegin() + d.size() - alignment_end));
    s.second.append(d.size() - alignment_end, '-');
  }

  //std::cout << "alignment_end, d.size() = " << alignment_end << "," << d.size() << "\n";

  auto add_del = [&]()
                 {
                   s.first.push_back(*std::next(d_begin, j - 1));
                   s.second.push_back('-');
                   --j;
                 };

  auto add_ins = [&]()
                 {
                   s.first.push_back('-');
                   s.second.push_back(*std::next(q_begin, i - 1));
                   --i;
                 };

  auto add_sub = [&]()
                 {
                   s.first.push_back(d[j - 1]);
                   s.second.push_back(q[i - 1]);
                   --i;
                   --j;
                 };

  while (i > 0 || j > 0)
  {
    std::size_t const v = j % t;
    std::size_t const e = j / t;

    if (i == 0)
    {
      while (j > 0)
        add_del();
    }
    else if (mB.is_del(i - 1, v, e))
    {
      while (j > 1 && mB.is_del_extend(i - 1, j % t, j / t))
        add_del();

      add_del();
    }
    else if (mB.is_ins(i - 1, v, e))
    {
      while (i > 1 && mB.is_ins_extend(i - 1, v, e))
        add_ins();

      add_ins();
    }
    else
    {
      add_sub();
    }
  }

  std::pair<std::string, std::string> out =
  {
    std::string(s.first.rbegin(), s.first.rend()),
    std::string(s.second.rbegin(), s.second.rend())
  };

  assert(out.first.size() == out.second.size());
  return out;
}


template <typename Tuint, typename Tit>
std::set<Event> inline
Aligner<Tuint, Tit>::get_edit_script(std::pair<std::string, std::string> const & s)
{
  using Tedit_script = std::set<Event>;

  Tedit_script edit_script;
  std::vector<char> s1;
  std::vector<char> s2;
  std::size_t pos = 0;

  auto are_sequences_empty = [&](){
                               return s1.size() == 0 && s2.size() == 0;
                             };

  auto add_to_edit_script = [&]()
                            {
                              Event new_edit =
                              {
                                pos - s1.size(),
                                std::string(s1.begin(), s1.end()),
                                std::string(s2.begin(), s2.end())
                              };

                              edit_script.insert(std::move(new_edit));

                              s1.clear();
                              s2.clear();
                            };


  for (std::size_t i = 0; i < s.first.size(); ++i)
  {
    // Both string cant have a gap
    assert(s.first[i] != '-' || s.second[i] != '-');

    if (s.first[i] == '-')
    {
      // Insertion
      if (are_sequences_empty() || (s2.size() > 0 && s1.size() == 0))
        s2.push_back(s.second[i]); // Extend the insertion
      else
        add_to_edit_script();
    }
    else if (s.second[i] == '-')
    {
      // Deletion
      if (are_sequences_empty() || (s1.size() > 0 && s2.size() == 0))
        s1.push_back(s.first[i]); // Extend the deletion
      else
        add_to_edit_script();

      ++pos;
    }
    else if (s.first[i] != s.second[i])
    {
      // Mismatch
      if (!are_sequences_empty())
        add_to_edit_script();

      ++pos;
      s1.push_back(s.first[i]);
      s2.push_back(s.second[i]);
      add_to_edit_script();
    }
    else
    {
      // Match
      if (!are_sequences_empty())
        add_to_edit_script();

      ++pos;
    }
  }

  if (!are_sequences_empty())
    add_to_edit_script();

  return edit_script;
}


template <typename Tuint, typename Tit>
std::vector<Cigar> inline
Aligner<Tuint, Tit>::get_cigar(std::pair<std::string, std::string> const & s)
{
  std::vector<Cigar> cigar;

  if (s.first.size() == 0)
    return cigar;

  auto get_op = [](char const c1, char const c2) -> CigarOperation
                {
                  if (c1 == '-')
                    return INSERTION;
                  else if (c2 == '-')
                    return DELETION;
                  else
                    return MATCH;
                };

  std::size_t prev_pos = 0;
  CigarOperation op = get_op(s.first[0], s.second[0]);

  for (std::size_t i = 1; i < s.first.size(); ++i)
  {
    // Both string cant have a gap
    assert(s.first[i] != '-' || s.second[i] != '-');
    CigarOperation const new_op = get_op(s.first[i], s.second[i]);

    if (op == new_op)
    {
      continue;
    }
    else
    {
      cigar.push_back({i - prev_pos, op});
      op = new_op;
      prev_pos = i;
    }
  }

  cigar.push_back({s.first.size() - prev_pos, op});
  return cigar;
}


} // namespace paw


#endif // IMPLEMENT_PAW
