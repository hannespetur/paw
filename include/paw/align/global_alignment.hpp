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
  long const right_v = m % t; // Vector that contains the rightmost element
  long const right_e = m / t; // The right-most element (in vector 'right_v')
  AlignmentResults<Tuint> ar;
  ar.mB = Backtrack<Tuint>(n, t);
  init_vH_up<Tuint>(opt.vH_up, opt.gap_open_val, std::numeric_limits<Tuint>::min());
  Tvec_pack vH(static_cast<std::size_t>(t), simdpp::make_int(2 * opt.gap_open_val + std::numeric_limits<Tuint>::min()));
  Tvec_pack vF(opt.vF_up);
  Tvec_pack vE(opt.vF_up);

#ifndef NDEBUG
  store_scores(opt, 0, vE);
#endif // NDEBUG

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

    if (opt.left_column_free)
    {
      Tarr_uint vF0;
      vF0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vF0[0], vF[0]);
      vF0[0] = simdpp::extract<0>(opt.vH_up[0]) + opt.y_gain;
      vF[0] = simdpp::load(&vF0[0]);
    }

    // In case right_v is 0
    if (opt.right_column_free && right_v == 0)
    {
      Tarr_uint vH_up_0;
      vH_up_0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vH_up_0[0], opt.vH_up[0]);

      Tarr_uint vF0;
      vF0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vF0[0], vF[0]);

      vF0[right_e] = vH_up_0[right_e] + opt.y_gain;
      vF[0] = simdpp::load(&vF0[0]);
    }

    ar.mB.set_ins_extend(i, 0, max_greater<Tuint>(vF[0], opt.vF_up[0]));
    ar.mB.set_ins(i, 0, max_greater<Tuint>(vH[0], vF[0]));
    /// Done calculating vector 0

    /// Calculate vectors v=1,...,t-1
    for (long v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH[v] = opt.vH_up[v - 1] + vW[v - 1];
      vF[v] = opt.vH_up[v] - gap_open_pack_y;

      // In case right_v is 0
      if (opt.right_column_free && right_v == v)
      {
        Tarr_uint vH_up_v;
        vH_up_v.fill(std::numeric_limits<Tuint>::min());
        simdpp::store_u(&vH_up_v[0], opt.vH_up[v]);

        Tarr_uint vF0;
        vF0.fill(std::numeric_limits<Tuint>::min());
        simdpp::store_u(&vF0[0], vF[v]);

        vF0[right_e] = vH_up_v[right_e] + opt.y_gain;
        vF[v] = simdpp::load(&vF0[0]);
      }

      ar.mB.set_ins_extend(i, v, max_greater<Tuint>(vF[v], opt.vF_up[v]));
      ar.mB.set_ins(i, v, max_greater<Tuint>(vH[v], vF[v]));
    } /// Done calculating vectors v=1,...,t-1

    {
      /// Deletions pass 1
      vE[0] = shift_one_right<Tuint>(vH[t - 1] - gap_open_pack_x, std::numeric_limits<Tuint>::min(), opt.reductions);

      for (long v = 1; v < t; ++v)
      {
        vE[v] = vH[v - 1] - gap_open_pack_x;
        ar.mB.set_del_extend(i, v, max_greater<Tuint>(vE[v], vE[v - 1]));
      }

      Tarr_uint vE0r;
      vE0r.fill(std::numeric_limits<Tuint>::min());
      Tarr_uint vE0;
      vE0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vE0[0], vE[0]);

      // Check for deletions in vector 0
      Tpack vE0r_pack = shift_one_right<Tuint>(vE[t - 1], std::numeric_limits<Tuint>::min(), opt.reductions);
      simdpp::store_u(&vE0r[0], vE0r_pack);
      bool is_improved = false;

      for (long e = 1; e < static_cast<long>(vE0r.size()); ++e)
      {
        long val = static_cast<long>(vE0r[e - 1]) + opt.reductions[e - 1] - opt.reductions[e];

        if (val > static_cast<long>(vE0r[e]))
          vE0r[e] = val;

        if (vE0r[e] > vE0[e])
        {
          vE0[e] = vE0r[e];
          is_improved = true;
        }
        else
        {
          vE0r[e] = vE0[e];
        }
      }

      if (is_improved)
        ar.mB.set_del_extend(i, 0, max_greater<Tuint>(vE[0], simdpp::load_u(&vE0r[0])));

      ar.mB.set_del(i, 0, max_greater<Tuint>(vH[0], vE[0]));
      /// Done with deletion pass 1

      /// Deletions pass 2
      if (is_improved)
      {
        for (long v = 1; v < t; ++v)
          ar.mB.set_del_extend(i, v, max_greater<Tuint>(vE[v], vE[v - 1]));
      }

      for (long v = 1; v < t; ++v)
        ar.mB.set_del(i, v, max_greater<Tuint>(vH[v], vE[v]));
      /// Done with deletion pass 2
    }

    std::swap(vF, opt.vF_up);
    std::swap(vH, opt.vH_up);

#ifndef NDEBUG
    store_scores(opt, i + 1l, vE);
#endif
  } /// End of outer loop

  Tvec_uint arr(S / sizeof(Tuint));
  simdpp::store_u(&arr[0], opt.vH_up[m % t]);
  ar.query_end = m;
  ar.database_end = n;

  ar.score = static_cast<long>(arr[m / t])
             + static_cast<long>(opt.reductions[m / t])
             - opt.top_left_score
             - n * opt.y_gain
             - m * opt.x_gain;
  ar.t = t;
  return ar;
}


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
