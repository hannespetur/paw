#pragma once

#include <cstdint>
#include <memory>
#include <type_traits>

#include <simdpp/simd.h>

#include <paw/align/alignment_cache.hpp>
#include <paw/align/alignment_results.hpp>
#include <paw/align/libsimdpp_utils.hpp>


namespace paw
{

template<typename Tuint>
struct AlignmentOptions
{
public:
  using Tpack = typename T<Tuint>::pack;

  bool left_column_free = false;
  bool right_column_free = false;
  bool continuous_alignment = false; /// When set, always continue with the same alignment as long as the query is the same

private:
  /// User options
  Tuint match = 2; /// Score of matches
  Tuint mismatch = 2; /// Penalty of mismatches
  Tuint gap_open = 5; /// Penalty of opening a gap
  Tuint gap_extend = 1; /// Penalty of extending a gap
  //bool is_traceback = true; /// Set if the alignment traceback is required

  // TODO: Implement usage of "convex" gap cost
  //Tuint gap_open_2 = 5; /// Penalty of opening a gap in when using a secondary value
  //Tuint gap_extend_2 = 1; /// Penalty of extending a gap in when using a secondary value
  ///

  /// Calculated values
  std::shared_ptr<AlignmentCache<Tuint> > ac{new AlignmentCache<Tuint>()}; /// Shared cache between alignments
  std::unique_ptr<AlignmentResults<Tuint> > ar{new AlignmentResults<Tuint>()}; /// Results of the alignment


public:

#ifndef NDEBUG
  // Store the score matrix in debug mode
  std::vector<std::vector<long> > score_matrix;
  std::vector<std::vector<long> > vE_scores;
  std::vector<std::vector<long> > vF_scores;
#endif // not NDEBUG

  AlignmentOptions() = default;


  /// \brief Sets score of a match
  /// It is assumed that the score is greater or equal to 0
  AlignmentOptions &
  set_match(long val)
  {
    match = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  /// \brief Sets penalty of mismatches
  /// It is assumed that the penalty is less or equal to 0
  AlignmentOptions &
  set_mismatch(long val)
  {
    mismatch = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  /// \brief Sets penalty of opening gaps
  /// It is assumed that the penalty is less or equal to 0
  AlignmentOptions &
  set_gap_open(long val)
  {
    gap_open = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  /// \brief Sets penalty of opening gaps
  /// It is assumed that the penalty is less or equal to 0
  AlignmentOptions &
  set_gap_extend(long val)
  {
    gap_extend = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  /*
  /// \brief Sets gap penalty for both opening and extending a gap
  /// It is assumed that the penalty is less or equal to 0
  AlignmentOptions &
  set_gap(long val)
  {
    set_gap_open(val);
    return set_gap_extend(val);
  }


  /// \brief Sets whether an alignment traceback is required
  AlignmentOptions &
  set_traceback(bool val)
  {
    is_traceback = val;
    return *this;
  }
  */


  inline Tuint
  get_gap_open_val_x() const
  {
    return gap_open - ac->x_gain;
  }


  inline Tuint
  get_gap_open_val_y() const
  {
    return gap_open - ac->y_gain;
  }


  inline Tuint get_match() const {return match;}
  inline Tuint get_mismatch() const {return mismatch;}
  inline Tuint get_gap_open() const {return gap_open;}
  inline Tuint get_gap_extend() const {return gap_extend;}
  inline AlignmentCache<Tuint> * get_alignment_cache() const {return ac.get();}
  inline AlignmentResults<Tuint> * get_alignment_results() const {return ar.get();}
  //inline bool get_is_traceback() const {return is_traceback;}

};


namespace SIMDPP_ARCH_NAMESPACE
{


template <typename Tuint>
inline void
init_vH_up(AlignmentOptions<Tuint> & opts)
{
  using Tvec_uint = typename T<Tuint>::vec_uint;

  long const gap_open_val = opts.get_alignment_cache()->gap_open_val;
  AlignmentResults<Tuint> & aln_results = *opts.get_alignment_results();
  Tvec_uint new_vH0(T<Tuint>::pack::length, 2 * gap_open_val + std::numeric_limits<Tuint>::min());
  assert(aln_results.vH_up.size() > 0);
  new_vH0[0] = gap_open_val * 3 + std::numeric_limits<Tuint>::min();
  aln_results.vH_up[0] = simdpp::load_u(&new_vH0[0]);
}


template<typename Tuint, typename Tseq>
void
set_query(AlignmentOptions<Tuint> & opt, Tseq const & seq)
{
  using Tpack = typename T<Tuint>::pack;
  using Tvec_pack = typename T<Tuint>::vec_pack;

  std::string new_query(begin(seq), end(seq));
  Tpack const min_value_pack = simdpp::make_int(std::numeric_limits<Tuint>::min());
  AlignmentCache<Tuint> & aln_cache = *opt.get_alignment_cache();
  AlignmentResults<Tuint> & aln_results = *opt.get_alignment_results();
  bool const is_new_query = new_query != aln_cache.query;

  if (is_new_query)
  {
    aln_cache.set_query(std::move(new_query));
    aln_cache.set_options(opt.get_match(),
                          opt.get_mismatch(),
                          opt.get_gap_open(),
                          opt.get_gap_extend());
  }

  if (!opt.continuous_alignment || is_new_query)
  {
    aln_results.vH_up = Tvec_pack(static_cast<std::size_t>(aln_cache.num_vectors),
                                  static_cast<Tpack>(simdpp::make_int(2 * aln_cache.gap_open_val + std::numeric_limits<Tuint>::min()))
      );

    init_vH_up(opt);
    aln_results.vF_up = Tvec_pack(static_cast<std::size_t>(aln_cache.num_vectors), min_value_pack);

    aln_results.reduction = - std::numeric_limits<Tuint>::min() - aln_cache.gap_open_val * 3;
    aln_results.reductions.fill(0);

#ifndef NDEBUG
    opt.score_matrix.clear();
    opt.vE_scores.clear();
    opt.vF_scores.clear();
#endif // NDEBUG
  }
}


template<typename Tuint>
void
reduce_too_high_scores(AlignmentOptions<Tuint> & opt)
{
  using Tpack = typename T<Tuint>::pack;
  using Tmask = typename T<Tuint>::mask;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  AlignmentCache<Tuint> const & aln_cache = *opt.get_alignment_cache();
  AlignmentResults<Tuint> & aln_results = *opt.get_alignment_results();
  long const num_vectors = aln_cache.num_vectors;

  if (simdpp::reduce_max(aln_results.vH_up[num_vectors - 1]) >= aln_cache.max_score_val)
  {
    /// Reducing value losslessly
    Tarr_uint vF0;
    vF0.fill(0);
    // Store the optimal scores in vector 0
    simdpp::store_u(&vF0[0], aln_results.vF_up[0]);
    assert(vF0.size() == aln_results.reductions.size());
    Tarr_uint new_reductions;
    new_reductions.fill(0);
    assert(vF0.size() == new_reductions.size());
    bool any_reductions = false;

    for (long e = 1; e < static_cast<long>(vF0.size()); ++e)
    {
      long new_reduction_val = static_cast<long>(vF0[e]) - static_cast<long>(2 * aln_cache.gap_open_val);

      if (new_reduction_val > 0)
      {
        aln_results.reductions[e] += new_reduction_val;
        new_reductions[e] = new_reduction_val;
        any_reductions = true;
      }
    }

    // Reduce values
    if (any_reductions)
    {
      Tpack new_reductions_pack = simdpp::load_u(&new_reductions[0]);

      for (long v = 0; v < num_vectors; ++v)
      {
        aln_results.vF_up[v] = aln_results.vF_up[v] - new_reductions_pack;
        aln_results.vH_up[v] = aln_results.vH_up[v] - new_reductions_pack;
      }
    }

    Tuint const max_score = simdpp::reduce_max(aln_results.vH_up[num_vectors - 1]);

    if (max_score >= aln_cache.max_score_val)
    {
      Tpack const two_gap_open_pack = simdpp::make_int(2 * aln_cache.gap_open_val);

      /// Reducing values lossily
      Tmask is_about_to_overflow = aln_results.vH_up[num_vectors - 1] >= static_cast<Tpack>(simdpp::make_int(aln_cache.max_score_val));

      Tpack const overflow_reduction =
        simdpp::blend(two_gap_open_pack,
                      static_cast<Tpack>(simdpp::make_zero()),
                      is_about_to_overflow
        );

      {
        Tarr_uint overflow_reduction_arr;
        overflow_reduction_arr.fill(0);
        simdpp::store_u(&overflow_reduction_arr[0], overflow_reduction);

        for (long e = 0; e < static_cast<long>(S / sizeof(Tuint)); ++e)
          aln_results.reductions[e] += overflow_reduction_arr[e];
      }

      for (long v = 0; v < num_vectors; ++v)
      {
        aln_results.vH_up[v] = simdpp::max(aln_results.vH_up[v] - overflow_reduction, two_gap_open_pack);
        aln_results.vF_up[v] = simdpp::max(aln_results.vF_up[v] - overflow_reduction, two_gap_open_pack);
      }
    }
  }
}


template<typename Tuint>
std::vector<long> inline
get_score_row(AlignmentOptions<Tuint> const & opt, long const i, typename T<Tuint>::vec_pack const & vX)
{
  using Tvec_uint = typename T<Tuint>::vec_uint;

  AlignmentCache<Tuint> const & aln_cache = *opt.get_alignment_cache();
  AlignmentResults<Tuint> const & aln_results = *opt.get_alignment_results();
  long const m = aln_cache.query_size;
  long const t = aln_cache.num_vectors;

  assert(t > 0l);
  assert(t == static_cast<long>(vX.size()));

  Tvec_uint vec(vX[0].length, 0);
  std::vector<Tvec_uint> mat(t, vec);
  std::vector<long> scores_row;
  scores_row.reserve(m + 1ul);

  for (long v = 0; v < t; ++v)
    simdpp::store_u(&mat[v][0], aln_results.vH_up[v]);

  for (long j = 0; j <= m; ++j)
  {
    long const v = j % t;
    long const e = j / t;
    assert(v < static_cast<long>(mat.size()));
    assert(e < static_cast<long>(mat[v].size()));

    long const adjustment = aln_results.reduction + aln_results.reductions[e] - aln_cache.y_gain * i - aln_cache.x_gain * j;
    scores_row.push_back(static_cast<long>(mat[v][e] + adjustment));
  }

  return scores_row;
}


template<typename Tuint>
AlignmentOptions<Tuint>
merge_score_matrices(std::vector<AlignmentOptions<Tuint>* > const & opts_vec_ptr)
{
  assert(opts_vec_ptr.size() > 2);

  AlignmentOptions<Tuint> final_opts(&(opts_vec_ptr.front()));
  AlignmentResults<Tuint> final_aln_results = *final_opts.get_alignment_results();
  std::vector<long> max_vH = get_score_row(final_opts, 0, final_aln_results.vH_up);
  std::vector<long> max_vF = get_score_row(final_opts, 0, final_aln_results.vF_up);

  auto set_max_scores = [&opts_vec_ptr](std::vector<long> & max_scores, std::vector<std::vector<long> > & merge_solver)
  {
    long const in_degree = static_cast<long>(opts_vec_ptr.size()); // InDegree of node that is getting merged
    long const num_scores = static_cast<long>(max_scores.size()); // Number of scores in a row
    merge_solver.clear();
    merge_solver.resize(max_scores.size(), std::vector<long>(1, 0));

    for (long i = 1; i < in_degree; ++i)
    {
      AlignmentOptions<Tuint> const & opt = &(opts_vec_ptr[i]);
      std::vector<long> scores = get_score_row(opt, 0, opt.get_alignment_results()->vH_up);
      assert (scores.size() == num_scores);

      for (long j = 0; j < num_scores; ++j)
      {
        if (scores[j] > max_scores[j])
        {
          // Score is higher than what was previously observed, mark is the best path
          max_scores[j] = scores[j];
          merge_solver[j].clear();
          merge_solver[j].push_back(i);
        }
        else if (scores[j] == max_scores[j])
        {
          // If score is equal then just add it among the best paths
          merge_solver[j].push_back(i);
        }
      }
    }
  };

  set_max_scores(max_vH, final_opts.merge_solver_vH);
  set_max_scores(max_vF, final_opts.merge_solver_vF);

  /// Initialize vH_up and vF_up such that minimum score is at least 3*opt.gap_open + min_value and the maximum score
  /// is at most opt.max_score_val
  {
    long const min_score_val = 3 * final_opts.gap_open + std::numeric_limits<Tuint>::min();
    long const max_score_val = final_opts.get_cache()->max_score_val;
    long const max_diff = max_score_val - min_score_val;

    std::array<long, S / sizeof(Tuint)> min_values;
    min_values.fill(std::numeric_limits<long>::max());

    std::array<long, S / sizeof(Tuint)> max_values;
    max_values.fill(std::numeric_limits<long>::min());

    assert(max_vH.size() - 1l == final_opts.get_cache()->query_size);
    assert(max_vF.size() - 1l == final_opts.get_cache()->query_size);

    long t = final_opts.get_cache()->num_vectors;

    for (long j = 0; j <= final_opts.get_cache()->query_size; ++j)
    {
      long const e = j / t;
      min_values[e] = std::min(min_values[e], std::min(max_vH[j], max_vF[j]));
      max_values[e] = std::max(max_values[e], std::max(max_vH[j], max_vF[j]));
    }
  }

  return final_opts;
}


#ifndef NDEBUG

template<typename Tuint>
inline void
store_scores(AlignmentOptions<Tuint> & opt,
             long const i,
             typename T<Tuint>::vec_pack const & vE
  )
{
  AlignmentResults<Tuint> const & aln_results = *opt.get_alignment_results();
  opt.score_matrix.push_back(get_score_row(opt, i, aln_results.vH_up));
  opt.vE_scores.push_back(get_score_row(opt, i, vE));
  opt.vF_scores.push_back(get_score_row(opt, i, aln_results.vF_up));
  assert(opt.score_matrix.size() == opt.score_matrix.size());
  assert(opt.vE_scores.size() == opt.vF_scores.size());
}

#endif // NDEBUG

} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
