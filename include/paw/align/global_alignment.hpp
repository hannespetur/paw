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
  using Tvec_pack = typename T<Tuint>::vec_pack;
  using Tvec_uint = typename T<Tuint>::vec_uint;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  paw::SIMDPP_ARCH_NAMESPACE::set_query<Tuint, Tseq>(opt, seq1);

  long const m = opt.query_size;
  long const t = opt.num_vectors; // Keep t as a local variable is it widely used
  long const n = std::distance(seq2.begin(), seq2.end());
  Tuint const top_left_score = opt.gap_open_val * 3 + std::numeric_limits<Tuint>::min();
  AlignmentResults<Tuint> ar;
  ar.mB = Backtrack<Tuint>(n, t);
  init_vH_up<Tuint>(opt.vH_up, opt.gap_open_val, std::numeric_limits<Tuint>::min());
  Tvec_pack vH(static_cast<std::size_t>(t), simdpp::make_int(2 * opt.gap_open_val + std::numeric_limits<Tuint>::min()));
  Tvec_pack vF(opt.vF_up);
  Tvec_pack vE(opt.vF_up);

#ifndef NDEBUG
  store_vH_up_scores(opt, m, 0);
#endif // NDEBUG

  Tpack const min_value_pack = simdpp::make_int(std::numeric_limits<Tuint>::min());
  Tpack const gap_open_pack_x = simdpp::make_int(opt.gap_open_val_x);
  Tpack const gap_open_pack_y = simdpp::make_int(opt.gap_open_val_y);

  /// Start of outer loop
  for (long i = 0; i < n; ++i)
  {
    reduce_too_high_scores(opt);

    // We need to increase fix vF_up if y_gain is more than gap_extend cost
    if (i > 0 && opt.y_gain > opt.get_gap_extend())
    {
      for (long v = 0; v < t; ++v)
        opt.vF_up[v] = opt.vF_up[v] + static_cast<Tpack>(simdpp::make_uint(opt.y_gain - opt.get_gap_extend()));
    }

    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Tvec_pack const & vW = opt.W_profile[magic_function(*std::next(seq2.begin(), i))];

    /// Calculate vector 0
    {
      auto const left = std::max(static_cast<Tuint>(simdpp::extract<0>(opt.vF_up[0])),
                                 static_cast<Tuint>(simdpp::extract<0>(opt.vH_up[0]) -
                                 opt.gap_open_val_y)
                                 );

      vH[0] = shift_one_right<Tuint>(opt.vH_up[t - 1] + vW[t - 1],
                                     left,
                                     opt.reductions);
    }

    // Check if any insertion have highest values
    vF[0] = opt.vH_up[0] - gap_open_pack_y;
    ar.mB.set_ins_extend(i, 0, max_greater<Tuint>(vF[0], opt.vF_up[0]));
    ar.mB.set_ins(i, 0, max_greater<Tuint>(vH[0], vF[0]));
    /// Done calculating vector 0

    vE[0] = min_value_pack;

    /// Calculate vectors v=1,...,t-1
    for (long v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH[v] = opt.vH_up[v - 1] + vW[v - 1];
      vF[v] = opt.vH_up[v] - gap_open_pack_y;
      ar.mB.set_ins_extend(i, v, max_greater<Tuint>(vF[v], opt.vF_up[v]));
      ar.mB.set_ins(i, v, max_greater<Tuint>(vH[v], vF[v]));

      /// Deletions pass 1
      vE[v] = vH[v - 1] - gap_open_pack_x;
      ar.mB.set_del_extend(i, v, max_greater<Tuint>(vE[v], vE[v - 1]));
    } /// Done calculating vectors v=1,...,t-1

    /// Deletions pass 2
    {
      // Calculate vE[0] after vH[t - 1] has been calculated
      Tpack const new_vE0 = vH[t - 1] - gap_open_pack_x;
      vE[0] = shift_one_right<Tuint>(new_vE0, std::numeric_limits<Tuint>::min(), opt.reductions);
      Tarr_uint vE0;
      vE0.fill(std::numeric_limits<Tuint>::min());

      /// Check for deletions in vector 0
      Tpack vE0_pack = shift_one_right<Tuint>(vE[t - 1], std::numeric_limits<Tuint>::min(), opt.reductions);
      simdpp::store_u(&vE0[0], vE0_pack);

      for (long e = 2; e < static_cast<long>(vE0.size()); ++e)
      {
        long val = static_cast<long>(vE0[e - 1]) + opt.reductions[e - 1] - opt.reductions[e];

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

    std::swap(vF, opt.vF_up);
    std::swap(vH, opt.vH_up);

#ifndef NDEBUG
    store_vH_up_scores(opt, m, i + 1l);
#endif
  } /// End of outer loop

  Tvec_uint arr(S / sizeof(Tuint));
  simdpp::store_u(&arr[0], opt.vH_up[m % t]);

  ar.score = static_cast<long>(arr[m / t])
             + static_cast<long>(opt.reductions[m / t])
             - top_left_score
             - n * opt.y_gain
             - m * opt.x_gain;

  //std::cout << "final score = " << static_cast<long>(arr[m / t]) << " + " << static_cast<long>(opt.reductions[m / t])
  //  << " - " << n * opt.y_gain << " - " << m * opt.x_gain << "\n";

  /*
  {
    long alignment_end = m;
    long i = n;
    long j = alignment_end;
    std::pair<std::string, std::string> s;

    // TODO: Fix this mess
    std::string d(q_begin, d_end);
    std::string q(d_begin, q_end);

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
      s.first.push_back(*std::next(q_begin, j - 1));
      s.second.push_back('-');
      --j;
      ++del_count;
    };

    auto add_del_ext = [&]()
    {
      //std::cerr << "DEL SELECTED" << j << "\n";
      s.first.push_back(*std::next(q_begin, j - 1));
      s.second.push_back('-');
      --j;
      ++del_ext_count;
    };

    auto add_ins = [&]()
    {
      //std::cerr << "INS SELECTED" << "\n";
      s.first.push_back('-');
      s.second.push_back(*std::next(d_begin, i - 1));
      --i;
      ++ins_count;
    };

    auto add_ins_ext = [&]()
    {
      //std::cerr << "INS SELECTED" << "\n";
      s.first.push_back('-');
      s.second.push_back(*std::next(d_begin, i - 1));
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
