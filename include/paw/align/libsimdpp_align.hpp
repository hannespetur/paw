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

std::vector<Cigar>
get_cigar(std::pair<std::string, std::string> const & s)
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


template <typename Tit, typename Tuint>
class Align
{
public:
  using Tpack = typename T<Tuint>::pack;
  using Tmask = typename T<Tuint>::mask;
  using Tvec_pack = typename T<Tuint>::vec_pack;
  using Tvec_uint = typename T<Tuint>::vec_uint;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  // p is the length (or cardinality) of the SIMD vectors
  std::size_t static const p = Tpack::length;


private:
  AlignerOptions<Tuint> const opt;
  Tit d_begin;
  Tit d_end;
  Tit q_begin;
  Tit q_end;

  long const m;
  long alignment_end = m;
  long n = 0;
  long const t;

  Backtrack<Tuint> mB; //(n /*n_row*/, t /*n_vectors in each score row*/);

  long calculate_scores();

public:
  Align(Tit _d_begin,
        Tit _d_end,
        AlignerOptions<Tuint> const & _opt = AlignerOptions<Tuint>(true /*default options*/)
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
    , mB()
  {}

  long align(Tit q_begin, Tit q_end);   // Align query
  std::pair<std::string, std::string> get_aligned_strings();
  //std::vector<Cigar> get_cigar(std::pair<std::string, std::string> const & s);

};


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
  //new_vH0.fill(gap_open_val);
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


template <typename Tit, typename Tuint>
inline void
calculate_DNA_W_profile(typename T<Tuint>::arr_vec_pack & W_profile,
                        Tit d_begin,
                        long const m,
                        Tuint const match,
                        Tuint const mismatch
                        )
{
  std::array<char, 4> constexpr DNA_BASES = {{'A', 'C', 'G', 'T'}};
  long const t = (m + T<Tuint>::pack::length) / T<Tuint>::pack::length;

  for (std::size_t i = 0; i < DNA_BASES.size(); ++i)
  {
    char const dna_base = DNA_BASES[i];
    auto & W = W_profile[i];
    W.reserve(t);

    for (long v = 0; v < t; ++v)
    {
      typename T<Tuint>::vec_uint seq(T<Tuint>::pack::length, mismatch);

      for (long e = 0, j = v; j < static_cast<long>(m); j += t, ++e)
      {
        if (dna_base == *(d_begin + j))
          seq[e] = match;
      }

      W.push_back(static_cast<typename T<Tuint>::pack>(simdpp::load_u(&seq[0])));
    }
  }

  assert(static_cast<std::size_t>(t) == W_profile[0].size());
  assert(static_cast<std::size_t>(t) == W_profile[1].size());
  assert(static_cast<std::size_t>(t) == W_profile[2].size());
  assert(static_cast<std::size_t>(t) == W_profile[3].size());
}


template <typename Tit, typename Tuint>
long
Align<Tit, Tuint>::align(Tit _q_begin, Tit _q_end)
{
  q_begin = _q_begin;
  q_end = _q_end;
  n = std::distance(q_begin, q_end);

  if (opt.backtracking)
  {
    //mB = Backtrack<Tuint>(n, t);
    if (n > static_cast<long>(mB.matrix.size()))
    {
      mB = Backtrack<Tuint>(n, t);
    }
    else
    {
      // Make sure the first n rows have only zeros
      for (long i = 0; i < n; ++i)
      {
        for (auto & v : mB.matrix[i])
          v = simdpp::make_zero();
      }
    }
  }

  return this->calculate_scores();
}


template <typename Tit, typename Tuint>
long
Align<Tit, Tuint>::calculate_scores()
{
  Tuint const x_gain = opt.gap_extend;
  Tuint const y_gain = std::max(static_cast<Tuint>(opt.gap_extend),
                                static_cast<Tuint>(opt.mismatch - x_gain));
  Tuint const gap_open_val_x = opt.gap_open - x_gain;
  Tuint const gap_open_val_y = opt.gap_open - y_gain;
  Tuint const gap_open_val = std::max(gap_open_val_x, gap_open_val_y);
  Tvec_pack vH_up(t, static_cast<Tpack>(simdpp::make_uint(gap_open_val)));
  init_vH_up<Tuint>(vH_up, gap_open_val);
  Tvec_pack vH(t, simdpp::make_uint(gap_open_val));
  Tvec_pack vF_up(t, simdpp::make_zero());
  Tvec_pack vF(vF_up);
  Tvec_pack vE(vF_up);
  Tuint const top_left_score = gap_open_val * 2;
  Tuint const max_score_val = std::numeric_limits<Tuint>::max() - 2 * gap_open_val;
  Tpack const max_score_pack = simdpp::make_uint(max_score_val);

  typename T<Tuint>::arr_vec_pack W_profile;

  calculate_DNA_W_profile<Tit, Tuint>(W_profile,
                                      d_begin,
                                      m,
                                      x_gain + y_gain + opt.match,
                                      x_gain + y_gain - opt.mismatch
                                      );


  std::array<long, S / sizeof(Tuint)> reductions;
  reductions.fill(0);
  Tpack gap_open_pack = simdpp::make_uint(gap_open_val);
  Tpack gap_open_pack_x = simdpp::make_uint(gap_open_val_x);

  /// Start of outer loop
  for (long i = 0; i < n; ++i)
  {
    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Tvec_pack const & vW = W_profile[magic_function(*std::next(q_begin, i))];

    /// Calculate vector 0
    {
      auto const left = std::max(static_cast<Tuint>(simdpp::extract<0>(vF_up[0])),
                                 static_cast<Tuint>(simdpp::extract<0>(vH_up[0]) -
                                                    gap_open_val));

      vH[0] = shift_one_right<Tuint>(vH_up[t - 1] + vW[t - 1],
                                     left,
                                     reductions);
    }

    // Check if any insertion have highest values
    vF[0] = vH_up[0] - gap_open_pack;
    mB.set_ins_extend(i, 0, max_greater<Tuint>(vF[0], vF_up[0]));
    mB.set_ins(i, 0, max_greater<Tuint>(vH[0], vF[0]));
    /// Done calculating vector 0

    //vE[0] = static_cast<Tpack>(simdpp::make_zero());

    /// Calculate vectors v=1,...,t-1
    for (long v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH[v] = vH_up[v - 1] + vW[v - 1];
      vF[v] = vH_up[v] - gap_open_pack;
      mB.set_ins_extend(i, v, max_greater<Tuint>(vF[v], vF_up[v]));
      mB.set_ins(i, v, max_greater<Tuint>(vH[v], vF[v]));

      // Deletions pass 1: Gap opens
      vE[v] = vH[v - 1] - gap_open_pack;
    } /// Done calculating vectors v=1,...,t-1

    /// Calculate vE[0] after vH[t - 1] has been calculated
    {
      Tpack const new_vE0 = vH[t - 1] - gap_open_pack_x;
      vE[0] = shift_one_right<Tuint>(new_vE0, 0, reductions);
    }

    /// Deletions pass 2: Gap extends
    {
      Tarr_uint vE0;
      vE0.fill(0);

      /// Two passes through the deletion vectors are required
      for (std::size_t c = 0; c < 2; ++c)
      {
        {
          Tpack vE0_pack = shift_one_right<Tuint>(vE[t - 1], 0, reductions);
          simdpp::store_u(&vE0[0], vE0_pack);
        }
        //simdpp::store_u(&vE0[0], static_cast<Tpack>(vE[0]));

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

          Tpack const vE0_pack = simdpp::load_u(&vE0[0]);
          mB.set_del_extend(i, 0, static_cast<Tmask>(max_greater<Tuint>(vE[0], vE0_pack)));
          mB.set_del(i, 0, static_cast<Tmask>(max_greater<Tuint>(vH[0], vE[0])));
        }

        /// Check for deletions in vectors 1,...,t-1
        for (std::size_t v = 1; v < vE.size(); ++v)
        {
          mB.set_del_extend(i, v, static_cast<Tmask>(max_greater<Tuint>(vE[v], vE[v - 1])));
        }
      }
    }

    // Check if any vE has higher scores than vH
    for (long v = 1; v < t; ++v)
      mB.set_del(i, v, max_greater<Tuint>(vH[v], vE[v]));

#ifndef NDEBUG
    //print_score_vectors(vH, vH_up, vE, vF, vF_up, vW); // Useful when debugging
#endif
    if (simdpp::reduce_max(vH[t - 1]) >= max_score_val)
    {
      Tarr_uint vF0;
      vF0.fill(0);
      // Store the optimal scores in vector 0
      simdpp::store_u(&vF0[0], vF[0]);
      assert(vF0.size() == reductions.size());
      Tarr_uint new_reductions;
      new_reductions.fill(0);
      assert(vF0.size() == new_reductions.size());
      bool any_reductions = false;

      for (long e = 1; e < static_cast<long>(vF0.size()); ++e)
      {
        //assert(vF0[e] + gap_open_val_x >= vF0[e - 1]);
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
        Tpack new_reductions_pack = simdpp::load_u(&new_reductions[0]);

        for (long v = 0; v < t; ++v)
        {
          vF[v] = vF[v] - new_reductions_pack;
          vH[v] = vH[v] - new_reductions_pack;
        }
      }

      Tpack const max_scores = vH[t - 1];
      Tuint const max_score = simdpp::reduce_max(max_scores);

      if (max_score >= max_score_val)
      {
        Tmask is_about_to_overflow = max_scores >= max_score_pack;
        Tpack overflow_reduction =
          simdpp::blend(gap_open_pack,
                        static_cast<Tpack>(simdpp::make_zero()),
                        is_about_to_overflow
                        );

        {
          Tarr_uint overflow_reduction_arr;
          overflow_reduction_arr.fill(0);
          simdpp::store_u(&overflow_reduction_arr[0], overflow_reduction);

          for (long e = 0; e < static_cast<long>(S / sizeof(Tuint)); ++e)
            reductions[e] += overflow_reduction_arr[e];
        }

        for (long v = 0; v < t; ++v)
        {
          vH[v] = simdpp::max(vH[v] - overflow_reduction, gap_open_pack);
          vF[v] = simdpp::max(vF[v] - overflow_reduction, gap_open_pack);
        }
      }
    }

    std::swap(vF, vF_up);
    std::swap(vH, vH_up);
  } /// End of outer loop

#ifndef NDEBUG
  //print_backtrack(mB);
#endif

  {
    alignment_end = m;
    Tvec_uint arr(S / sizeof(Tuint));
    simdpp::store_u(&arr[0], vH_up[m % t]);
    return arr[m / t]
           + reductions[m / t]
           - top_left_score
           - n * y_gain
           - m * x_gain;
  }
}


template <typename Tit, typename Tuint>
std::pair<std::string, std::string>
Align<Tit, Tuint>::get_aligned_strings()
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


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
