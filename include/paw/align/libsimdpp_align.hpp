#pragma once

#include <array> // std::array
#include <cassert>
#include <cstdint> // uint8_t
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <paw/align/aligner_options.hpp>
#include <paw/align/event.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>

#include <simdpp/simd.h>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tit>
class Align
{
public:
  // p is the length (or cardinality) of the SIMD vectors
  std::size_t static const p = T::pack::length;
  void calculate_DNA_W_profile();


private:
  AlignerOptions<T::uint> const opt;
  //std::vector<Event> free_snp_edits;
  //std::vector<Event> free_del_edits;
  //std::vector<Event> free_ins_edits;
  Tit d_begin;
  Tit d_end;
  Tit q_begin;
  Tit q_end;

  long const m = 0;
  long alignment_end = m;
  long n = 0;
  long const t = 0;

  T::uint const x_gain;
  T::uint const y_gain;

  T::uint const gap_open_val_x;
  T::uint const gap_open_val_y;
  T::uint const gap_open_val;

  T::uint const max_score_val;

  std::array<long, S / sizeof(T::uint)> reductions;
  long total_reductions = 0;
  Row vH_up; // Previous H row
  Row vH;    // Current H row
  Row vE;    // Current E row
  Row vF_up; // Previous F row
  Row vF;    // Current F row
  T::arr_row W_profile;
  Backtrack mB; //(n /*n_row*/, t /*n_vectors in each score row*/);
  T::uint top_left_score = 0;

  void calculate_scores();
  void check_gap_extend_deletions();
  void check_gap_extend_deletions_with_backtracking(std::size_t const i,
                                                    std::array<long, S / sizeof(T::uint)> const &
                                                    );
  void init_score_vectors();

public:
  Align(Tit _d_begin,
        Tit _d_end,
        AlignerOptions<T::uint> const & _opt = AlignerOptions<T::uint>(true /*default options*/),
        std::set<Event> const & /*free_edits*/ = std::set<Event>()
        )
    : opt(_opt)
    , d_begin(_d_begin)
    , d_end(_d_end)
    , q_begin(_d_begin)
    , q_end(_d_end)
    , m(std::distance(_d_begin, _d_end))
    , alignment_end(m)
    , n(0)
    , t((m + p) / p)
    , x_gain(opt.gap_extend)
    , y_gain(std::max(opt.gap_extend, static_cast<T::uint>(_opt.mismatch - x_gain)))
    , gap_open_val_x(opt.gap_open - x_gain)
    , gap_open_val_y(opt.gap_open - y_gain)
    , gap_open_val(std::max(gap_open_val_x, gap_open_val_y))
    , max_score_val(std::numeric_limits<T::uint>::max() - 2 * gap_open_val)
    , reductions({})
    , total_reductions(0)
    , vH_up(m + 1, gap_open_val)
    , vH(m + 1, gap_open_val)
    , vE(m + 1)
    , vF_up(m + 1)
    , vF(m + 1)
    , W_profile
    {
      {
        Row(m + 1),
        Row(m + 1),
        Row(m + 1),
        Row(m + 1)
      }
    }
    , mB()
    , top_left_score(0)
  {
    /*
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
    */
  }

  int64_t align(Tit q_begin, Tit q_end);   // Align query
  std::pair<std::string, std::string> get_aligned_strings();
  std::set<Event> get_edit_script(std::pair<std::string, std::string> const & s);
  std::vector<Cigar> get_cigar(std::pair<std::string, std::string> const & s);

  inline T::arr_row const & get_W_profile()
  {
    return W_profile;
  }
};


inline long
magic_function(char const c)
{
  return 0x03 & ((c >> 2) ^ (c >> 1));
}


inline T::mask
max_greater(T::pack & v1, T::pack const & v2)
{
  // TODO: Only use simdpp::max when no backtracking.
  T::mask is_greater = v2 > v1;
  v1 = simdpp::blend(v2, v1, is_greater);
  return is_greater;
}


template <typename Tit>
void
Align<Tit>::calculate_DNA_W_profile()
{
  std::array<char, 4> constexpr DNA_BASES = {{'A', 'C', 'G', 'T'}};
  T::uint const match = x_gain + y_gain + opt.match;
  T::uint const mismatch = x_gain + y_gain - opt.mismatch;
  long const t = W_profile[0].vectors.size();
  assert(t > 0);

  for (std::size_t i = 0; i < DNA_BASES.size(); ++i)
  {
    char const dna_base = DNA_BASES[i];
    auto & W = W_profile[i];

    for (long v = 0; v < t; ++v)
    {
      T::vec_uint seq(p, mismatch);

      for (long e = 0, j = v; j < static_cast<long>(m); j += t, ++e)
      {
        if (dna_base == *(d_begin + j))
          seq[e] = match;
      }

      W.vectors[v] = static_cast<T::pack>(simdpp::load_u(&seq[0]));
    }

    // Update the W_profile if there are any free mismatches
    /*
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
    */
  }
}


template <typename Tit>
void
Align<Tit>::init_score_vectors()
{
  // For less verbosity of code
  auto const x = gap_open_val;

  //if (opt.top_row_gap_open_free)
  //{
  //  for (auto & v : vH_up.vectors)
  //    v = simdpp::make_int(x * 2);
  //
  //  for (auto & v : vH.vectors)
  //    v = simdpp::make_int(x * 2);
  //}
  //else
  {
    //for (auto & v : vH_up.vectors)
    //  v = simdpp::make_int(x);
    T::arr_uint new_vH0;
    new_vH0.fill(x);
    new_vH0[0] = x * 2;
    vH_up.vectors[0] = simdpp::load_u(&new_vH0[0]);
    //simdpp::insert<0>(vH_up.vectors[0], x * 2);
    //vH_up.vectors[0] = simdpp::make_uint(x * 2, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x);

    //for (auto & v : vH.vectors)
    //  v = simdpp::make_int(x);
  }

  //for (auto & v : vF_up.vectors)
  //  v = simdpp::make_zero();

  // TODO: Implement
  //if (opt.top_row_free)
  //{
  //  for (std::size_t v = 0; v < t; ++v)
  //  {
  //    for (std::size_t e = 0, j = v; j < m; ++e, j += t)
  //    {
  //      Tuint const gain = j * x_gain;
  //      vH_up.vectors[v][e] = 2 * gap_open_val + gain;
  //      vH.vectors[v][e] = 2 * gap_open_val + gain;
  //      vF_up.vectors[v][e] = gap_open_val + gain;
  //    }
  //  }
  //}

  top_left_score = x * 2;
}


template <typename Tit>
int64_t
Align<Tit>::align(Tit _q_begin, Tit _q_end)
{
  q_begin = _q_begin;
  q_end = _q_end;
  n = std::distance(q_begin, q_end);
  this->init_score_vectors();
  total_reductions = 0;

  if (opt.backtracking)
  {
    mB = Backtrack(n, t);
    //if (n > mB.matrix.size())
    //{
    //  mB = Backtrack<(n, t);
    //}
    //else
    //{
    //  // Make sure the first n rows have only zeros
    //  for (std::size_t i = 0; i < n; ++i)
    //  {
    //    for (auto & v : mB.matrix[i])
    //      v = simdpp::make_zero();
    //  }
    //}
  }

  this->calculate_scores();

  /*
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
  */
  {
    alignment_end = m;
    T::vec_uint arr(S / sizeof(T::uint));
    simdpp::store_u(&arr[0], vH_up.vectors[m % t]);
    return arr[m / t]
           + reductions[m / t]
           - top_left_score
           - n * y_gain
           - m * x_gain;
  }
}


template <typename Tit>
void
Align<Tit>::calculate_scores()
{
  reductions.fill(0);
  //T::pack reductions_diff_pack = simdpp::make_zero();
  T::pack gap_open_pack = simdpp::make_uint(gap_open_val);
  T::pack gap_open_pack_x = simdpp::make_uint(gap_open_val_x);
  /*T::uint const max_gain_per_row = std::max(static_cast<T::uint>(0),
                                            static_cast<T::uint>(opt.match + x_gain + y_gain)
                                            );

  // Keep track of current max score
  T::uint current_max_score = top_left_score; //simdpp::extract<0>(vH_up.vectors[0]);
  */

  /// Start of outer loop
  for (long i = 0; i < n; ++i)
  {
    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Row const & vW = W_profile[magic_function(*std::next(q_begin, i))];

    // Clear vE vectors
    for (auto & vE_vec : vE.vectors)
      vE_vec = simdpp::make_zero();

    /// Calculate vector 0
    {
      auto const left = std::max(static_cast<T::uint>(simdpp::extract<0>(vF_up.vectors[0])),
                                 static_cast<T::uint>(simdpp::extract<0>(vH_up.vectors[0]) - gap_open_val));

      vH.vectors[0] = shift_one_right(vH_up.vectors[t - 1] + vW.vectors[t - 1], left, reductions);
    }

    // Check if any insertion have highest values
    //if (backtracking)
    {
      vF.vectors[0] = vH_up.vectors[0] - gap_open_pack;

      //if (left_column_gap_open_free)
      //  vF.vectors[0][0] = vH_up.vectors[0][0];
      //
      //if (right_column_gap_open_free && right_v == 0)
      //  vF.vectors[0][right_e] = vH_up.vectors[0][right_e];

      mB.set_ins_extend(i, 0, max_greater(vF.vectors[0], vF_up.vectors[0]));
      mB.set_ins(i, 0, max_greater(vH.vectors[0], vF.vectors[0]));
    }
    //else
    //{
    //  vF.vectors[0] = vH_up.vectors[0] - gap_open_pack;
    //
    //  if (left_column_gap_open_free)
    //    vF.vectors[0][0] = vH_up.vectors[0][0];
    //
    //  if (right_column_gap_open_free && right_v == 0)
    //    vF.vectors[0][right_e] = vH_up.vectors[0][right_e];
    //
    //  vF.vectors[0] = boost::simd::max(vF.vectors[0], vF_up.vectors[0]);
    //  vH.vectors[0] = boost::simd::max(vH.vectors[0], vF.vectors[0]);
    //}
    /// Done calculating vector 0

    /// Calculate vectors 1,...,v-1
    for (long v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH.vectors[v] = vH_up.vectors[v - 1] + vW.vectors[v - 1];

      //if (backtracking)
      {
        vF.vectors[v] = vH_up.vectors[v] - gap_open_pack;

        /// Slight performance gain if this if statement is moved out of the loop
        //if (right_column_gap_open_free && v == right_v)
        //  vF.vectors[right_v][right_e] = vH_up.vectors[right_v][right_e];

        mB.set_ins_extend(i, v, max_greater(vF.vectors[v], vF_up.vectors[v]));
        mB.set_ins(i, v, max_greater(vH.vectors[v], vF.vectors[v]));
      }
      //else
      //{
      //  vF.vectors[v] = boost::simd::max(vH_up.vectors[v] - gap_open_pack, vF_up.vectors[v]);
      //
      //  // Slight performance gain if this if statement is moved out of the loop
      //  if (right_column_gap_open_free && v == right_v)
      //    vF.vectors[right_v][right_e] = vH_up.vectors[right_v][right_e];
      //
      //  vH.vectors[v] = boost::simd::max(vH.vectors[v], vF.vectors[v]);
      //}

      // Deletions pass 1: Gap opens
      vE.vectors[v] = vH.vectors[v - 1] - gap_open_pack;
    } /// Done calculating vectors 1,...,v-1


    /// Calculate vE.vector[0]
    {
      T::pack const new_vE0 = vH.vectors[t - 1] - gap_open_pack_x;
      vE.vectors[0] = shift_one_right(new_vE0, 0, reductions);
    }

    //vE.vectors[0] = simdpp::move8_r<1>(vH.vectors[t - 1] - gap_open_pack_x);

    /// Deletions pass 2: Gap extends
    //if (backtracking)
    {
      check_gap_extend_deletions_with_backtracking(i, reductions);

      // Check if any vE has higher scores than vH
      for (long v = 0; v < t; ++v)
        mB.set_del(i, v, max_greater(vH.vectors[v], vE.vectors[v]));
    }
    //else
    //{
    //  check_gap_extend_deletions();
    //
    //  // Check if any vE has higher scores than vH
    //  for (std::size_t v = 0; v < t; ++v)
    //    vH.vectors[v] = boost::simd::max(vH.vectors[v], vE.vectors[v]);
    //}

#ifndef NDEBUG
    print_score_vectors(vH, vH_up, vE, vF, vF_up, vW); // Useful when debugging
#endif
    //std::cout << "max score = " << static_cast<uint64_t>(current_max_score) << "\n";

    //if (i % 100 == 99)
    {
      T::arr_uint vF0;
      vF0.fill(0);
      // Store the optimal scores in vector 0
      simdpp::store_u(&vF0[0], vF.vectors[0]);
      assert(vF0.size() == reductions.size());
      T::arr_uint new_reductions;
      new_reductions.fill(0);
      assert(vF0.size() == new_reductions.size());
      bool any_reductions = false;

      for (long e = 1; e < static_cast<long>(vF0.size()); ++e)
      {
        assert(vF0[e] + gap_open_val_x >= vF0[e - 1]);
        long new_reduction_val = static_cast<long>(vF0[e]) - static_cast<long>(2 * gap_open_val_x);

        if (new_reduction_val > 0)
        {
          new_reductions[e] = new_reduction_val;
          any_reductions = true;
          reductions[e] += static_cast<long>(new_reductions[e]);
        }
      }

      // Reduce values
      if (any_reductions)
      {
        T::pack new_reductions_pack = simdpp::load_u(&new_reductions[0]);

        for (long v = 0; v < t; ++v)
        {
          vF.vectors[v] = vF.vectors[v] - new_reductions_pack;
          vH.vectors[v] = vH.vectors[v] - new_reductions_pack;
        }
      }

      T::pack const max_scores = vH.vectors[t - 1];
      T::uint const max_score = simdpp::reduce_max(max_scores);

      if (max_score > max_score_val)
      {
        T::pack max_score_pack = simdpp::make_uint(max_score_val);
        T::mask is_about_to_overflow = max_scores > max_score_pack;
        T::pack overflow_reduction = simdpp::blend(gap_open_pack, static_cast<T::pack>(simdpp::make_zero()), is_about_to_overflow);

        {
          T::arr_uint overflow_reduction_arr;
          overflow_reduction_arr.fill(0);
          simdpp::store_u(&overflow_reduction_arr[0], overflow_reduction);

          for (long e = 0; e < static_cast<long>(S / sizeof(T::uint)); ++e)
            reductions[e] += overflow_reduction_arr[e];
        }

        for (long v = 0; v < t; ++v)
        {
          vH.vectors[v] = simdpp::max(vH.vectors[v] - overflow_reduction, gap_open_pack);
          vF.vectors[v] = simdpp::max(vF.vectors[v] - overflow_reduction, gap_open_pack);
        }
      }
    }

    std::swap(vF.vectors, vF_up.vectors);
    std::swap(vH.vectors, vH_up.vectors);
  } /// End of outer loop

#ifndef NDEBUG
  //print_backtrack(mB);
#endif
}


template <typename Tit>
void
Align<Tit>::check_gap_extend_deletions_with_backtracking(std::size_t const i, std::array<long, S / sizeof(T::uint)> const & reductions)
{
  T::arr_uint vE0;
  vE0.fill(0);

  /// Two passes through the deletion vectors are required
  for (std::size_t c = 0; c < 2; ++c)
  {
    {
      T::pack vE0_pack = shift_one_right(vE.vectors[t - 1], 0, reductions);
      simdpp::store_u(&vE0[0], vE0_pack);
    }

    /// Check for deletions in vector 0
    {
      for (long e = 2; e < static_cast<long>(vE0.size()); ++e)
      {
        long val = static_cast<long>(vE0[e - 1]) + reductions[e - 1] - reductions[e];

        if (val > 0)
        {
          vE0[e] = std::max(val, static_cast<long>(vE0[e]));
        }
      }

      T::pack const vE0_pack = simdpp::load_u(&vE0[0]);
      T::mask const del_extend_mask_0 = max_greater(vE.vectors[0], vE0_pack);
      mB.set_del_extend(i, 0, del_extend_mask_0);
    }

    /// Check for deletions in vectors 1,...,t-1
    for (std::size_t v = 1; v < vE.vectors.size(); ++v)
    {
      T::mask const del_extend_mask_v = max_greater(vE.vectors[v], vE.vectors[v - 1]);
      mB.set_del_extend(i, v, del_extend_mask_v);
    }
  }
}


template <typename Tit>
std::pair<std::string, std::string>
Align<Tit>::get_aligned_strings()
{
  long i = n;
  long j = alignment_end;
  std::pair<std::string, std::string> s;

  // TODO: Fix this mess
  std::string d(d_begin, d_end);
  std::string q(q_begin, q_end);

  long del_count = 0;
  long ins_count = 0;
  long del_ext_count = 0;
  long ins_ext_count = 0;
  long match_count = 0;
  long mismatch_count = 0;

  if (alignment_end < static_cast<long>(d.size()))
  {
    s.first.append(std::string(d.rbegin(), d.rbegin() + d.size() - alignment_end));
    s.second.append(d.size() - alignment_end, '-');

    del_count = 1;
    del_ext_count = d.size() - alignment_end - 1;
  }

  //std::cout << "alignment_end, d.size() = " << alignment_end << "," << d.size() << "\n";

  auto add_del = [&]()
                 {
                   //std::cerr << "DEL SELECTED" << j << "\n";
                   s.first.push_back(*std::next(d_begin, j - 1));
                   s.second.push_back('-');
                   --j;
                   ++del_count;
                 };

  auto add_del_ext = [&]()
                 {
                   //std::cerr << "DEL SELECTED" << j << "\n";
                   s.first.push_back(*std::next(d_begin, j - 1));
                   s.second.push_back('-');
                   --j;
                   ++del_ext_count;
                 };

  auto add_ins = [&]()
                 {
                   //std::cerr << "INS SELECTED" << "\n";
                   s.first.push_back('-');
                   s.second.push_back(*std::next(q_begin, i - 1));
                   --i;
                   ++ins_count;
                 };

  auto add_ins_ext = [&]()
                 {
                   //std::cerr << "INS SELECTED" << "\n";
                   s.first.push_back('-');
                   s.second.push_back(*std::next(q_begin, i - 1));
                   --i;
                   ++ins_ext_count;
                 };

  auto add_sub = [&]()
                 {
                   //std::cerr << "SUB SELECTED" << "\n";
                   if (d[j - 1] == q[i - 1])
                   {
                     ++match_count;
                   }
                   else
                   {
                     ++mismatch_count;
                   }

                   s.first.push_back(d[j - 1]);
                   s.second.push_back(q[i - 1]);
                   --i;
                   --j;


                 };

  while (i > 0 || j > 0)
  {
    if (j == 0)
    {
      while (i > 1)
        add_ins_ext();

      add_ins();
      break;
    }

    assert(i >= 0);
    assert(j >= 0);
    //std::cerr << "i,j = " << i << "," << j << std::endl;
    std::size_t const v = j % t;
    std::size_t const e = j / t;

    if (i == 0)
    {
      assert(j > 0);
      add_del();

      while (j > 0)
        add_del_ext();
    }
    else if (mB.is_del(i - 1, v, e))
    {
      while (j > 1 && mB.is_del_extend(i - 1, j % t, j / t))
        add_del_ext();

      assert(j > 0);
      add_del();
    }
    else if (mB.is_ins(i - 1, v, e))
    {
      while (i > 1 && mB.is_ins_extend(i - 1, v, e))
        add_ins_ext();

      assert(i > 0);
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


template <typename Tit>
std::vector<Cigar>
Align<Tit>::get_cigar(std::pair<std::string, std::string> const & s)
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


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
