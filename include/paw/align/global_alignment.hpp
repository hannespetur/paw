#pragma once

#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
#include <paw/align/event.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>


#include <simdpp/simd.h>

#include <array> // std::array
#include <cassert> // assert
#include <cstdint> // uint8_t, ...
#include <iostream> // std::cerr
#include <limits> // std::numeric_limits
#include <numeric>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tseq, typename Tuint>
AlignmentResults<Tuint>
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opt
                 )
{
  using Tpack = typename T<Tuint>::pack;
  using Tmask = typename T<Tuint>::mask;
  using Tvec_pack = typename T<Tuint>::vec_pack;
  using Tvec_uint = typename T<Tuint>::vec_uint;
  using Tarr_uint = typename T<Tuint>::arr_uint;
  using Tarr_vec_pack = typename T<Tuint>::arr_vec_pack;

  auto d_begin = begin(seq1);
  auto d_end = end(seq1);
  auto q_begin = begin(seq2);
  auto q_end = end(seq2);

  long const m = std::distance(d_begin, d_end);
  long const t = (m + Tpack::length) / Tpack::length;
  long const n = std::distance(q_begin, q_end);
  assert(t >= 0);
  AlignmentResults<Tuint> ar;
  ar.mB = Backtrack<Tuint>(n, t);
  Tuint const x_gain = opt.get_gap_extend();
  Tuint const y_gain = std::max(static_cast<Tuint>(opt.get_gap_extend()),
                                static_cast<Tuint>(opt.get_mismatch() - x_gain));
  Tuint const gap_open_val_x = opt.get_gap_open() - x_gain;
  Tuint const gap_open_val_y = opt.get_gap_open() - y_gain;
  Tuint const gap_open_val = std::max(gap_open_val_x, gap_open_val_y);
  Tuint const min_value = std::numeric_limits<Tuint>::min();
  Tpack const min_value_pack = simdpp::make_int(min_value);
  Tvec_pack vH_up(static_cast<std::size_t>(t), static_cast<Tpack>(simdpp::make_int(2 * gap_open_val + min_value)));
  init_vH_up<Tuint>(vH_up, gap_open_val, min_value);
  Tvec_pack vH(static_cast<std::size_t>(t), simdpp::make_int(2 * gap_open_val + min_value));
  Tvec_pack vF_up(static_cast<std::size_t>(t), min_value_pack);
  Tvec_pack vF(vF_up);
  Tvec_pack vE(vF_up);
  Tuint const match = x_gain + y_gain + opt.get_match();
  Tuint const mismatch = x_gain + y_gain - opt.get_mismatch();
  Tuint const max_score_val = std::numeric_limits<Tuint>::max() - match - gap_open_val;
  Tpack const max_score_pack = simdpp::make_int(max_score_val);
  std::array<long, S / sizeof(Tuint)> reductions;
  reductions.fill(0);
  Tpack const two_gap_open_pack = simdpp::make_int(2 * gap_open_val);
  Tpack const gap_open_pack_x = simdpp::make_int(gap_open_val_x);
  Tpack const gap_open_pack_y = simdpp::make_int(gap_open_val_y);
  Tuint const top_left_score = gap_open_val * 3 + min_value;

  Tarr_vec_pack W_profile;
#ifndef NDEBUG
  //std::cout << "vH TOP"; print_score_vector_standard<Tuint>(m, vH_up, top_left_score, x_gain, 0, reductions);
#endif // NDEBUG

  /// Calculate DNA W_profile
  {
    std::array<char, 4> constexpr DNA_BASES = {{'A', 'C', 'G', 'T'}};

    for (std::size_t i = 0; i < DNA_BASES.size(); ++i)
    {
      char const dna_base = DNA_BASES[i];
      auto & W = W_profile[i];
      W.reserve(t);

      for (long v = 0; v < t; ++v)
      {
        Tvec_uint seq(T<Tuint>::pack::length, mismatch);

        for (long e = 0, j = v; j < m; j += t, ++e)
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
  } /// Done calculating DNA W_profile

  /// Start of outer loop
  for (long i = 0; i < n; ++i)
  {
    // We need to increase fix vF_up if y_gain is more than gap_extend cost
    if (i > 0 && y_gain > opt.get_gap_extend())
    {
      for (long v = 0; v < t; ++v)
        vF_up[v] = vF_up[v] + static_cast<Tpack>(simdpp::make_uint(y_gain - opt.get_gap_extend()));
    }

    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Tvec_pack const & vW = W_profile[magic_function(*std::next(q_begin, i))];

    /// Calculate vector 0
    {
      auto const left = std::max(static_cast<Tuint>(simdpp::extract<0>(vF_up[0])),
                                 static_cast<Tuint>(simdpp::extract<0>(vH_up[0]) -
                                                    gap_open_val_y));

      vH[0] = shift_one_right<Tuint>(vH_up[t - 1] + vW[t - 1],
                                     left,
                                     reductions);
    }

    // Check if any insertion have highest values
    vF[0] = vH_up[0] - gap_open_pack_y;
    ar.mB.set_ins_extend(i, 0, max_greater<Tuint>(vF[0], vF_up[0]));
    ar.mB.set_ins(i, 0, max_greater<Tuint>(vH[0], vF[0]));
    /// Done calculating vector 0

    vE[0] = min_value_pack;

    /// Calculate vectors v=1,...,t-1
    for (long v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH[v] = vH_up[v - 1] + vW[v - 1];
      vF[v] = vH_up[v] - gap_open_pack_y;
      ar.mB.set_ins_extend(i, v, max_greater<Tuint>(vF[v], vF_up[v]));
      ar.mB.set_ins(i, v, max_greater<Tuint>(vH[v], vF[v]));

      /// Deletions pass 1
      vE[v] = vH[v - 1] - gap_open_pack_x;
      ar.mB.set_del_extend(i, v, max_greater<Tuint>(vE[v], vE[v - 1]));
    } /// Done calculating vectors v=1,...,t-1

    /// Deletions pass 2
    {
      // Calculate vE[0] after vH[t - 1] has been calculated
      Tpack const new_vE0 = vH[t - 1] - gap_open_pack_x;
      vE[0] = shift_one_right<Tuint>(new_vE0, min_value, reductions);
      Tarr_uint vE0;
      vE0.fill(min_value);

      /// Check for deletions in vector 0
      Tpack vE0_pack = shift_one_right<Tuint>(vE[t - 1], min_value, reductions);
      simdpp::store_u(&vE0[0], vE0_pack);

      for (long e = 2; e < static_cast<long>(vE0.size()); ++e)
      {
        long val = static_cast<long>(vE0[e - 1]) + reductions[e - 1] - reductions[e];

        if (val > static_cast<long>(vE0[e]))
          vE0[e] = val;
      }

      {
        Tpack const new_vE0_pack = simdpp::load_u(&vE0[0]);
        ar.mB.set_del_extend(i, 0, max_greater<Tuint>(vE[0], new_vE0_pack));
        ar.mB.set_del(i, 0, max_greater<Tuint>(vH[0], vE[0]));

        /// Check for deletions in vectors 1,...,t-1
        for (long v = 1; v < t; ++v)
        {
          ar.mB.set_del_extend(i, v, max_greater<Tuint>(vE[v], vE[v - 1]));
          ar.mB.set_del(i, v, max_greater<Tuint>(vH[v], vE[v]));
        }
      }
    }

    if (simdpp::reduce_max(vH[t - 1]) >= max_score_val)
    {
      /// Reducing value losslessly
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
        long new_reduction_val = static_cast<long>(vF0[e]) - static_cast<long>(2 * gap_open_val);

        if (new_reduction_val > 0)
        {
          reductions[e] += new_reduction_val;
          new_reductions[e] = new_reduction_val;
          any_reductions = true;
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
        /// Reducing values lossily
        Tmask is_about_to_overflow = max_scores >= max_score_pack;
        Tpack overflow_reduction =
          simdpp::blend(two_gap_open_pack,
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
          vH[v] = simdpp::max(vH[v] - overflow_reduction, two_gap_open_pack);
          vF[v] = simdpp::max(vF[v] - overflow_reduction, two_gap_open_pack);
        }
      }
    }

#ifndef NDEBUG
    //std::cout << "vE FIN"; print_score_vector_standard<Tuint>(m, vE, top_left_score, x_gain, (i + 1l) * y_gain, reductions);
    //std::cout << "vH    "; print_score_vector_standard<Tuint>(m, vH, top_left_score, x_gain, (i + 1l) * y_gain, reductions);

    //print_score_vectors<Tuint>(m, vH, vH_up, vE, vF, vF_up, vW, top_left_score, x_gain, i * y_gain, reductions); // Useful when debugging
#endif

    std::swap(vF, vF_up);
    std::swap(vH, vH_up);
  } /// End of outer loop

  Tvec_uint arr(S / sizeof(Tuint));
  simdpp::store_u(&arr[0], vH_up[m % t]);

  ar.score = static_cast<long>(arr[m / t])
             + static_cast<long>(reductions[m / t])
             - static_cast<long>(top_left_score)
             - n * y_gain
             - m * x_gain;

  ar.last_vF = std::move(vF_up);
  ar.last_vH = std::move(vH_up);


  /*
  {
    long alignment_end = m;
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

    if (m < static_cast<long>(d.size()))
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
      else if (ar.mB.is_del(i - 1, v, e))
      {
        while (j > 1 && ar.mB.is_del_extend(i - 1, j % t, j / t))
          add_del_ext();

        assert(j > 0);
        add_del();
      }
      else if (ar.mB.is_ins(i - 1, v, e))
      {
        while (i > 1 && ar.mB.is_ins_extend(i - 1, v, e))
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

    std::cout << "\n" << out.first << "\n" << out.second << std::endl;
  }
  */

  return ar;
}


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
