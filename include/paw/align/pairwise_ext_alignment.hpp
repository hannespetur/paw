#pragma once

#include <paw/align/alignment_ext_cache.hpp>
#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
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

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
                   Tseq const & seq2,
                   AlignmentOptions<Tuint> & opts);


/// See namespaces at:
/// https://github.com/p12tic/libsimdpp/blob/3ab0be0e5aa0773f152d7d759400173d64253534/simdpp/detail/insn_id.h

namespace arch_null
{

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
                   Tseq const & seq2,
                   AlignmentOptions<Tuint> & opts);

}

namespace arch_sse2
{

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
                   Tseq const & seq2,
                   AlignmentOptions<Tuint> & opts);

}

namespace arch_sse3
{

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
                   Tseq const & seq2,
                   AlignmentOptions<Tuint> & opts);

}

namespace arch_sse4p1
{

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
                   Tseq const & seq2,
                   AlignmentOptions<Tuint> & opts);

}

namespace arch_sse4p1_popcnt
{

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
                   Tseq const & seq2,
                   AlignmentOptions<Tuint> & opts);

}

namespace arch_popcnt_avx
{

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
                   Tseq const & seq2,
                   AlignmentOptions<Tuint> & opts);

}

namespace arch_popcnt_avx2
{

template <typename Tseq, typename Tuint>
void
pairwise_ext_alignment(Tseq const & seq1,
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
pairwise_ext_alignment(Tseq const & seq1, // seq1 is query
                       Tseq const & seq2, // seq2 is database
                       AlignmentOptions<Tuint> & opt)
{
  using Tpack = typename T<Tuint>::pack;
  using Tvec_pack = typename T<Tuint>::vec_pack;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  AlignmentExtCache<Tuint> aln_cache; // The ExtCache forces x_gain == 0
  paw::SIMDPP_ARCH_NAMESPACE::set_query_ext<Tuint, Tseq>(opt, aln_cache, seq1);
  paw::SIMDPP_ARCH_NAMESPACE::set_database_ext<Tuint, Tseq>(aln_cache, seq2);

  assert(opt.get_alignment_results());
  AlignmentResults & aln_results = *opt.get_alignment_results();

  int const m = aln_cache.query_size; // Local variable for the query size
  int const t = aln_cache.num_vectors; // Keep t as a local variable as it is widely used
  int const n = std::distance(seq2.begin(), seq2.end());
  int max_score{std::numeric_limits<int>::min()};
  int max_score_i{-1};
  int max_score_v{-1};
  int max_score_e{-1};
  int const right_v = m % t; // Vector that contains the rightmost element
  int const right_e = m / t; // The right-most element (in vector 'right_v')
  Tpack const gap_open_pack_x = simdpp::make_int(opt.get_gap_open());
  Tuint const gap_open_val_y = opt.get_gap_open_val_y(aln_cache); // Store once
  Tpack const gap_open_pack_y = simdpp::make_int(gap_open_val_y);

  Tvec_pack vH;
  vH.reserve(t);
  //(static_cast<std::size_t>(t), simdpp::make_int(2 * aln_cache.gap_open_val));

  for (long v{0}; v < t; ++v)
  {
    vH.push_back(aln_cache.vH_up[v] - gap_open_pack_x);
  }

  Tvec_pack vF(aln_cache.vF_up);
  Tvec_pack vE(aln_cache.vF_up);

#ifndef NDEBUG
  store_scores(opt, aln_cache, 0, vE);
#endif // NDEBUG

  /// Start of outer loop
  for (long i{0}; i < n; ++i)
  {
    //aln_cache.reduce_every_element(aln_cache.y_gain);
    // reduce_too_high_scores(aln_cache); // TODO make a similar function that must reduce every element

    // We need to fix vF_up if y_gain is more than gap_extend cost
    if (i > 0 && aln_cache.y_gain > opt.get_gap_extend())
    {
      for (long v{0}; v < t; ++v)
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

      vH[0] = shift_one_right<Tuint>(aln_cache.vH_up[t - 1] + vW[t - 1], left);
    }

    // Check if any insertion have highest values
    vF[0] = aln_cache.vH_up[0] - gap_open_pack_y;

    aln_cache.mB.set_ins_extend(i, 0, max_greater<Tuint>(vF[0], aln_cache.vF_up[0]));
    aln_cache.mB.set_ins(i, 0, max_greater<Tuint>(vH[0], vF[0]));
    /// Done calculating vector 0

    /// Calculate vectors v=1,...,t-1
    for (long v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH[v] = aln_cache.vH_up[v - 1] + vW[v - 1];
      vF[v] = aln_cache.vH_up[v] - gap_open_pack_y;

      aln_cache.mB.set_ins_extend(i, v, max_greater<Tuint>(vF[v], aln_cache.vF_up[v]));
      aln_cache.mB.set_ins(i, v, max_greater<Tuint>(vH[v], vF[v]));
    } /// Done calculating vectors v=1,...,t-1

    {
      /// Deletions pass 1
      vE[0] = shift_one_right<Tuint>(vH[t - 1] - gap_open_pack_x, std::numeric_limits<Tuint>::min());

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
      Tpack vE0r_pack = shift_one_right<Tuint>(vE[t - 1], std::numeric_limits<Tuint>::min());
      simdpp::store_u(&vE0r[0], vE0r_pack);
      bool is_improved = false;

      for (long e = 1; e < static_cast<long>(vE0r.size()); ++e)
      {
        long val = static_cast<long>(vE0r[e - 1]);

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

    // get the max scores
    if ((i + 1) < n)
    {
      int current_max_vector{right_v};

      // reduce by clip
      int current_max_score = static_cast<int>(simdpp::reduce_max(aln_cache.vH_up[right_v] + aln_cache.vX[right_v])) -
        static_cast<int>(opt.get_clip());

      int score_right;

      // first check the last query element
      {
        std::vector<Tuint> arr_right(S / sizeof(Tuint));
        simdpp::store_u(&arr_right[0], aln_cache.vH_up[right_v] + aln_cache.vX[right_v]);
        score_right = arr_right[right_e];
        //std::cerr << " score_right=" << score_right << '\n';

        if (score_right >= current_max_score)
        {
          current_max_score = score_right; // no clip reduction
          /* current_max_vector = right_v; */ // commented because it is already set
        }
      }

      for (long v{0}; v < t; ++v)
      {
        if (v == right_v)
          continue;

        auto score = static_cast<int>(simdpp::reduce_max(aln_cache.vH_up[v] + aln_cache.vX[v])) - static_cast<int>(opt.get_clip());

        if (score > current_max_score)
        {
          current_max_score = score;
          current_max_vector = v;
        }
      }

      int const reduction = (i + 1) * aln_cache.y_gain - aln_cache.reduction;
      int const corrected_max_score = current_max_score - reduction;

      if (corrected_max_score >= max_score)
      {
        // get which element in the current_max_vector had the maximum score
        // std::cerr << " new max corrected score=" << corrected_max_score << " from score=" << current_max_score << '\n';
        std::vector<Tuint> arr(S / sizeof(Tuint));
        simdpp::store_u(&arr[0], aln_cache.vH_up[current_max_vector] + aln_cache.vX[current_max_vector]);
        auto find_max_it = std::find(arr.begin(), arr.end(), static_cast<Tuint>(current_max_score + opt.get_clip()));

        if (find_max_it == arr.end())
        {
          // happens only if right_ve is best alignment
          assert(current_max_vector == right_v);
          assert(current_max_score == score_right);

          max_score_v = right_v; // set which vector contains the maximum score
          max_score_e = right_e; // set which element has the maximum score
        }
        else
        {
          max_score_v = current_max_vector; // set which vector contains the maximum score
          max_score_e = std::distance(arr.begin(), find_max_it); // set which element has the maximum score
        }

        max_score = corrected_max_score; // set a new max_score to the current one
        max_score_i = i; // set which row contains the maximum score
      }
      else
      {
        // check if it is impossible to get a better score
        int const potential_score = corrected_max_score + opt.get_match() * (n - 1 - i) + opt.get_clip();

        if (potential_score < max_score)
        {
          // std::cerr << " potential new score too low=" << potential_score << " < " << max_score << '\n';
          // std::cerr << " max_score_i,v,e=" << max_score_i << ", " << max_score_v << ", " << max_score_e << '\n';
          // std::cerr << " t=" << t << '\n';
          break;
        }
      }

      // check if the score is too high
      if ((current_max_score + opt.get_clip()) >= aln_cache.max_score_val)
      {
        // std::cerr << "Triggered a reduction of values" << std::endl;
        // reduce all scores by gap_open_val + match_val
        Tpack const reduce_pack = simdpp::make_int(aln_cache.gap_open_val + aln_cache.match_val);
        Tpack const reduce_pack_2x = simdpp::make_int(2 * (aln_cache.gap_open_val + aln_cache.match_val));
        aln_cache.reduction += aln_cache.gap_open_val + aln_cache.match_val;

        for (int v{0}; v < t; ++v)
        {
          aln_cache.vH_up[v] = simdpp::max(aln_cache.vH_up[v], reduce_pack_2x) - reduce_pack;
          aln_cache.vF_up[v] = simdpp::max(aln_cache.vH_up[v], reduce_pack_2x) - reduce_pack;
        }
      }
    }
  } /// End of outer loop

  std::vector<Tuint> arr(S / sizeof(Tuint));
  simdpp::store_u(&arr[0], aln_cache.vH_up[m % t] + aln_cache.vX[m % t]);

  /*
#ifndef NDEBUG
  if (aln_cache.num_vectors == 1)
  {
    std::cerr << "last index: " << (m/t) << '\n';
    std::cerr << "row: ";

    for (int e{0}; e <= aln_cache.query_size; ++e)
    {
      std::cerr << static_cast<int>(arr[e]) << ",";
    }

    std::cerr << '\n';
  }
#endif // NDEBUG
  */

  aln_results.query_end = m;
  aln_results.database_end = n;
  aln_results.score = static_cast<long>(arr[m / t])
    + static_cast<long>(aln_cache.reduction)
//                    - m * static_cast<long>(aln_cache.x_gain)
                    - n * static_cast<long>(aln_cache.y_gain);

  /*
#ifndef NDEBUG
  std::cerr << static_cast<int>(arr[m/t]) << " + " //
            << aln_cache.reduction << " - " //
//            << (m * static_cast<long>(aln_cache.x_gain)) << " - "
            << (n * static_cast<long>(aln_cache.y_gain)) << " = "
            << aln_results.score << '\n';
#endif // NDEBUG
  */

  if (max_score > static_cast<int>(aln_results.score))
  {
    //std::cerr << "clip at: " << max_score_i + 1
    //          << " with score=" << max_score << " > " << aln_results.score << std::endl;
    aln_results.query_end = t * max_score_e + max_score_v;
    aln_results.database_end = max_score_i + 1;
    aln_results.score = max_score;
  }

  aln_results.clip_begin = aln_results.query_begin;
  aln_results.clip_end = aln_results.query_end;

  if (opt.get_aligned_strings)
    aln_results.get_aligned_strings(aln_cache, seq1, seq2);

  if (opt.get_cigar_string)
    aln_results.get_cigar_string(aln_cache);
}


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
