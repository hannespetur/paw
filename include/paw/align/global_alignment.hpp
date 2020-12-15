#pragma once

#include <paw/align/alignment_cache.hpp>
#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
//#include <paw/align/event.hpp>
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

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);


/// See namespaces at:
/// https://github.com/p12tic/libsimdpp/blob/3ab0be0e5aa0773f152d7d759400173d64253534/simdpp/detail/insn_id.h

namespace arch_null
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);

}

namespace arch_sse2
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);

}

namespace arch_sse3
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);

}

namespace arch_sse4p1
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);

}

namespace arch_sse4p1_popcnt
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);

}

namespace arch_popcnt_avx
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);

}

namespace arch_popcnt_avx2
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts);

}

} // namespace paw


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1, // seq1 is query
                 Tseq const & seq2, // seq2 is database
                 AlignmentOptions<Tuint> & opt)
{
  using Tpack = typename T<Tuint>::pack;
  using Tvec_pack = typename T<Tuint>::vec_pack;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  AlignmentCache<Tuint> aln_cache;
  paw::SIMDPP_ARCH_NAMESPACE::set_query<Tuint, Tseq>(opt, aln_cache, seq1);
  paw::SIMDPP_ARCH_NAMESPACE::set_database<Tuint, Tseq>(aln_cache, seq2);

  assert(opt.get_alignment_results());
  AlignmentResults<Tuint> & aln_results = *opt.get_alignment_results();

  long const m = aln_cache.query_size; // Local variable for the query size
  long const t = aln_cache.num_vectors; // Keep t as a local variable is it widely used
  long const right_v = m % t; // Vector that contains the rightmost element
  long const right_e = m / t; // The right-most element (in vector 'right_v')
  long const n = std::distance(seq2.begin(), seq2.end());

  Tvec_pack vH(static_cast<std::size_t>(t), simdpp::make_int(
                 2 * aln_cache.gap_open_val + std::numeric_limits<Tuint>::min()));
  Tvec_pack vF(aln_cache.vF_up);
  Tvec_pack vE(aln_cache.vF_up);

#ifndef NDEBUG
  store_scores(opt, aln_cache, 0, vE);
#endif // NDEBUG

  Tpack const gap_open_pack_x = simdpp::make_int(opt.get_gap_open_val_x(aln_cache));
  Tuint const gap_open_val_y = opt.get_gap_open_val_y(aln_cache); // Store once
  Tpack const gap_open_pack_y = simdpp::make_int(gap_open_val_y);


  /// Start of outer loop
  for (long i = 0; i < n; ++i)
  {
    reduce_too_high_scores(aln_cache);

    // We need to increase fix vF_up if y_gain is more than gap_extend cost
    if (i > 0 && aln_cache.y_gain > opt.get_gap_extend())
    {
      for (long v = 0; v < t; ++v)
      {
        aln_cache.vF_up[v] = aln_cache.vF_up[v] +
                             static_cast<Tpack>(simdpp::make_uint(aln_cache.y_gain -
                                                                  opt.get_gap_extend()));
      }
    }

    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Tvec_pack const & vW = aln_cache.W_profile[magic_function(*std::next(seq2.begin(), i))];

    /// Calculate vector 0
    {
      auto const left = std::max(static_cast<Tuint>(simdpp::extract<0>(aln_cache.vF_up[0])),
                                 static_cast<Tuint>(simdpp::extract<0>(aln_cache.vH_up[0]) -
                                                    gap_open_val_y)
                                 );

      vH[0] = shift_one_right<Tuint>(aln_cache.vH_up[t - 1] + vW[t - 1],
                                     left,
                                     aln_cache.reductions);
    }

    // Check if any insertion have highest values
    vF[0] = aln_cache.vH_up[0] - gap_open_pack_y;

    if (opt.left_column_free)
    {
      Tarr_uint vF0;
      vF0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vF0[0], vF[0]);
      vF0[0] = simdpp::extract<0>(aln_cache.vH_up[0]) + aln_cache.y_gain;
      vF[0] = simdpp::load(&vF0[0]);
    }

    // In case right_v is 0
    if (opt.right_column_free && right_v == 0)
    {
      Tarr_uint vH_up_0;
      vH_up_0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vH_up_0[0], aln_cache.vH_up[0]);

      Tarr_uint vF0;
      vF0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vF0[0], vF[0]);

      vF0[right_e] = vH_up_0[right_e] + aln_cache.y_gain;
      vF[0] = simdpp::load(&vF0[0]);
    }

    aln_cache.mB.set_ins_extend(i, 0, max_greater<Tuint>(vF[0], aln_cache.vF_up[0]));
    aln_cache.mB.set_ins(i, 0, max_greater<Tuint>(vH[0], vF[0]));
    /// Done calculating vector 0

    /// Calculate vectors v=1,...,t-1
    for (long v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH[v] = aln_cache.vH_up[v - 1] + vW[v - 1];
      vF[v] = aln_cache.vH_up[v] - gap_open_pack_y;

      // In case right_v is 0
      if (opt.right_column_free && right_v == v)
      {
        Tarr_uint vH_up_v;
        vH_up_v.fill(std::numeric_limits<Tuint>::min());
        simdpp::store_u(&vH_up_v[0], aln_cache.vH_up[v]);

        Tarr_uint vF0;
        vF0.fill(std::numeric_limits<Tuint>::min());
        simdpp::store_u(&vF0[0], vF[v]);

        vF0[right_e] = vH_up_v[right_e] + aln_cache.y_gain;
        vF[v] = simdpp::load(&vF0[0]);
      }

      aln_cache.mB.set_ins_extend(i, v, max_greater<Tuint>(vF[v], aln_cache.vF_up[v]));
      aln_cache.mB.set_ins(i, v, max_greater<Tuint>(vH[v], vF[v]));
    } /// Done calculating vectors v=1,...,t-1

    {
      /// Deletions pass 1
      vE[0] = shift_one_right<Tuint>(vH[t - 1] - gap_open_pack_x,
                                     std::numeric_limits<Tuint>::min(), aln_cache.reductions);

      for (long v = 1; v < t; ++v)
      {
        vE[v] = vH[v - 1] - gap_open_pack_x;
        aln_cache.mB.set_del_extend(i, v, max_greater<Tuint>(vE[v], vE[v - 1]));
      }

      Tarr_uint vE0r;
      vE0r.fill(std::numeric_limits<Tuint>::min());
      Tarr_uint vE0;
      vE0.fill(std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vE0[0], vE[0]);

      // Check for deletions in vector 0
      Tpack vE0r_pack = shift_one_right<Tuint>(vE[t - 1],
                                               std::numeric_limits<Tuint>::min(),
                                               aln_cache.reductions);
      simdpp::store_u(&vE0r[0], vE0r_pack);
      bool is_improved = false;

      for (long e = 1; e < static_cast<long>(vE0r.size()); ++e)
      {
        long val = static_cast<long>(vE0r[e - 1]) + aln_cache.reductions[e - 1] -
                   aln_cache.reductions[e];

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
        aln_cache.mB.set_del_extend(i, 0, max_greater<Tuint>(vE[0], simdpp::load_u(&vE0r[0])));

      aln_cache.mB.set_del(i, 0, max_greater<Tuint>(vH[0], vE[0]));
      /// Done with deletion pass 1

      /// Deletions pass 2
      if (is_improved)
      {
        for (long v = 1; v < t; ++v)
          aln_cache.mB.set_del_extend(i, v, max_greater<Tuint>(vE[v], vE[v - 1]));
      }

      for (long v = 1; v < t; ++v)
        aln_cache.mB.set_del(i, v, max_greater<Tuint>(vH[v], vE[v]));
      /// Done with deletion pass 2
    }

    std::swap(vF, aln_cache.vF_up);
    std::swap(vH, aln_cache.vH_up);

#ifndef NDEBUG
    store_scores(opt, aln_cache, i + 1l, vE);
#endif
  } /// End of outer loop

  std::vector<Tuint> arr(S / sizeof(Tuint));
  simdpp::store_u(&arr[0], aln_cache.vH_up[m % t]);
  aln_results.query_end = m;
  aln_results.database_end = n;
  aln_cache.reduce_every_element(-n * aln_cache.y_gain);
  aln_results.score = static_cast<long>(arr[m / t])
                      + static_cast<long>(aln_cache.reductions[m / t])
                      - m * aln_cache.x_gain;

  if (opt.get_aligned_strings)
  {
    aln_results.get_aligned_strings(aln_cache, seq1, seq2);
  }
}


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
