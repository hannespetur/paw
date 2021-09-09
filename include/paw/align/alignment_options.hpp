#pragma once

#include <cstdint>
#include <memory>
#include <set>
#include <type_traits>

#include <simdpp/simd.h>

#include <paw/align/alignment_cache.hpp>
#include <paw/align/alignment_results.hpp>
//#include <paw/align/event.hpp>
#include <paw/align/libsimdpp_utils.hpp>
#include <paw/align/sequence_utils.hpp>


namespace paw
{

class Event2;

template <typename Tuint>
struct AlignmentOptions
{
public:
  bool left_column_free{false};
  bool right_column_free{false};
  bool continuous_alignment{false}; /// When set, always continue with the same alignment as long as the query is the same
  bool get_aligned_strings{false};
  bool get_cigar_string{false};
  bool is_clip{false};
  bool is_query_reverse_complement{false};
  std::set<Event2> free_edits; // free SNP events

private:
  /// User options
  Tuint match = 2; /// Score of matches
  Tuint mismatch = 2; /// Penalty of mismatches
  Tuint gap_open = 5; /// Penalty of opening a gap
  Tuint gap_extend = 1; /// Penalty of extending a gap
  Tuint clip = 5; /// Penalty of clipping the query
  //bool is_traceback = true; /// Set if the alignment traceback is required

  // TODO: Implement usage of "convex" gap cost
  //Tuint gap_open_2 = 5; /// Penalty of opening a gap in when using a secondary value
  //Tuint gap_extend_2 = 1; /// Penalty of extending a gap in when using a secondary value
  ///

  /// Calculated values
  //std::shared_ptr<AlignmentCache<Tuint> > ac{new AlignmentCache<Tuint>()}; /// Shared cache between alignments
  std::unique_ptr<AlignmentResults<Tuint> > ar{new AlignmentResults<Tuint>()}; /// Results of the alignment


public:

#ifndef NDEBUG
  /// DEBUG MODE ONLY: Store the score matrix
  std::vector<std::vector<long> > score_matrix;
  std::vector<std::vector<long> > vE_scores;
  std::vector<std::vector<long> > vF_scores;
#endif // NDEBUG


  AlignmentOptions()
    : left_column_free(false)
    , right_column_free(false)
    , continuous_alignment(false)
    , get_aligned_strings(false)
    , get_cigar_string(false)
    , is_clip(false)
    , is_query_reverse_complement(false)
    , match(2)
    , mismatch(2)
    , gap_open(5)
    , gap_extend(1)
    , clip(5)
    , ar(new AlignmentResults<Tuint>())
  {}


  AlignmentOptions(AlignmentOptions const & ao)
    : left_column_free(ao.left_column_free)
    , right_column_free(ao.right_column_free)
    , continuous_alignment(ao.continuous_alignment)
    , match(ao.match)
    , mismatch(ao.mismatch)
    , gap_open(ao.gap_open)
    , gap_extend(ao.gap_extend)
    , clip(ao.clip)
    , ar(new AlignmentResults<Tuint>())
  {
  }


  AlignmentOptions(AlignmentOptions && ao) noexcept
    : left_column_free(ao.left_column_free)
    , right_column_free(ao.right_column_free)
    , continuous_alignment(ao.continuous_alignment)
    , match(ao.match)
    , mismatch(ao.mismatch)
    , gap_open(ao.gap_open)
    , gap_extend(ao.gap_extend)
    , clip(ao.clip)
    , ar(std::forward<std::unique_ptr<AlignmentResults<Tuint> > >(ao.ar))
  {
    //ar = std::move(ao.ar);
  }


  AlignmentOptions<Tuint> &
  operator=(AlignmentOptions const & ao)
  {
    left_column_free = ao.left_column_free;
    right_column_free = ao.right_column_free;
    continuous_alignment = ao.continuous_alignment;

    match = ao.match;
    mismatch = ao.mismatch;
    gap_open = ao.gap_open;
    gap_extend = ao.gap_extend;
    clip = ao.clip;

    ar = std::unique_ptr<AlignmentResults<Tuint> >(new AlignmentResults<Tuint>());
    return *this;
  }


  AlignmentOptions<Tuint> &
  operator=(AlignmentOptions && ao) noexcept
  {
    left_column_free = ao.left_column_free;
    right_column_free = ao.right_column_free;
    continuous_alignment = ao.continuous_alignment;

    match = ao.match;
    mismatch = ao.mismatch;
    gap_open = ao.gap_open;
    gap_extend = ao.gap_extend;
    clip = ao.clip;

    ar = std::move(ao.ar);
    return *this;
  }


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


  /// \brief Sets penalty of clipping
  /// It is assumed that the penalty is less or equal to 0
  AlignmentOptions &
  set_clip(long val)
  {
    clip = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
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
  get_gap_open_val_x(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> const & aln_cache) const
  {
    return gap_open - aln_cache.x_gain;
  }


  inline Tuint
  get_gap_open_val_y(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> const & aln_cache) const
  {
    return gap_open - aln_cache.y_gain;
  }


  inline Tuint
  get_match() const {return match;}
  inline Tuint
  get_mismatch() const {return mismatch;}
  inline Tuint
  get_gap_open() const {return gap_open;}
  inline Tuint
  get_gap_extend() const {return gap_extend;}
  inline Tuint
  get_clip() const {return clip;}

  //inline AlignmentCache<Tuint> *
  //get_alignment_cache() const {return ac.get();}

  inline AlignmentResults<Tuint> *
  get_alignment_results() const {return ar.get();}
  //inline bool get_is_traceback() const {return is_traceback;}

};


namespace SIMDPP_ARCH_NAMESPACE
{


template <typename Tuint, typename Tseq>
void
set_query(AlignmentOptions<Tuint> & opt, AlignmentCache<Tuint> & aln_cache, Tseq const & seq)
{
  using Tpack = typename T<Tuint>::pack;
  using Tvec_pack = typename T<Tuint>::vec_pack;

  std::string new_query;

  if (!opt.is_query_reverse_complement)
  {
    new_query = std::string(seq.data(), seq.size()); // TODO this allocation is unneccesary
  }
  else
  {
    new_query.resize(seq.size(), '\0');
    std::transform(seq.rbegin(), seq.rend(), new_query.begin(), paw::complement);
  }

  Tpack const min_value_pack = simdpp::make_int(std::numeric_limits<Tuint>::min());
  //bool const is_new_query = new_query != aln_cache.query;

  //if (is_new_query)
  {
    aln_cache.set_query(std::move(new_query));
    aln_cache.set_options(opt.get_match(),
                          opt.get_mismatch(),
                          opt.get_gap_open(),
                          opt.get_gap_extend());
  }

  // set free snps
  for (Event2 const & e : opt.free_edits)
  {
    assert(e.is_snp());
    aln_cache.set_free_snp(e.pos, e.alt[0]);
  }

  aln_cache.vH_up = Tvec_pack(static_cast<std::size_t>(aln_cache.num_vectors),
                              static_cast<Tpack>(simdpp::make_int(2 * aln_cache.gap_open_val +
                                                                  std::numeric_limits<Tuint>::min()))
                              );

  // init vH up
  {
    long const gap_open_val = aln_cache.gap_open_val;
    std::vector<Tuint> new_vH0(T<Tuint>::pack::length,
                               2 * gap_open_val + std::numeric_limits<Tuint>::min());

    assert(aln_cache.vH_up.size() > 0);
    new_vH0[0] = gap_open_val * 3 + std::numeric_limits<Tuint>::min();
    aln_cache.vH_up[0] = simdpp::load_u(&new_vH0[0]);

  }

  aln_cache.vF_up = Tvec_pack(static_cast<std::size_t>(aln_cache.num_vectors), min_value_pack);
  aln_cache.reductions.fill(static_cast<long>(-std::numeric_limits<Tuint>::min()) - aln_cache.gap_open_val * 3);


#ifndef NDEBUG
  opt.score_matrix.clear();
  opt.vE_scores.clear();
  opt.vF_scores.clear();
#endif // NDEBUG
}


template <typename Tuint, typename Tseq>
void
set_database(AlignmentCache<Tuint> & aln_cache, Tseq const & seq)
{
  aln_cache.mB = Backtrack<Tuint>(std::distance(begin(seq), end(seq)), aln_cache.num_vectors);
}


template <typename Tuint>
void
reduce_too_high_scores(AlignmentCache<Tuint> & aln_cache)
{
  using Tpack = typename T<Tuint>::pack;
  using Tmask = typename T<Tuint>::mask;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  //AlignmentCache<Tuint> const & aln_cache = *opt.get_alignment_cache();
  //AlignmentResults<Tuint> & aln_results = *opt.get_alignment_results();
  long const num_vectors = aln_cache.num_vectors;

  if (simdpp::reduce_max(aln_cache.vH_up[num_vectors - 1]) >= aln_cache.max_score_val)
  {
    /// Reducing value losslessly
    Tarr_uint vF0;
    vF0.fill(0);
    // Store the optimal scores in vector 0
    simdpp::store_u(&vF0[0], aln_cache.vF_up[0]);
    assert(vF0.size() == aln_cache.reductions.size());
    Tarr_uint new_reductions;
    new_reductions.fill(0);
    assert(vF0.size() == new_reductions.size());
    bool any_reductions = false;

    for (long e = 1; e < static_cast<long>(vF0.size()); ++e)
    {
      long new_reduction_val = static_cast<long>(vF0[e]) - static_cast<long>(2 * aln_cache.gap_open_val);

      if (new_reduction_val > 0)
      {
        aln_cache.reductions[e] += new_reduction_val;
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
        aln_cache.vF_up[v] = aln_cache.vF_up[v] - new_reductions_pack;
        aln_cache.vH_up[v] = aln_cache.vH_up[v] - new_reductions_pack;
      }
    }

    Tuint const max_score = simdpp::reduce_max(aln_cache.vH_up[num_vectors - 1]);

    if (max_score >= aln_cache.max_score_val)
    {
      Tpack const two_gap_open_pack = simdpp::make_int(2 * aln_cache.gap_open_val);

      /// Reducing values lossily
      Tmask is_about_to_overflow = aln_cache.vH_up[num_vectors - 1] >=
                                   static_cast<Tpack>(simdpp::make_int(aln_cache.max_score_val));

      Tpack const overflow_reduction =
        simdpp::blend(two_gap_open_pack,
                      static_cast<Tpack>(simdpp::make_zero()),
                      is_about_to_overflow
                      );

      {
        Tarr_uint overflow_reduction_arr;
        overflow_reduction_arr.fill(0);
        simdpp::store_u(&overflow_reduction_arr[0], overflow_reduction);

        for (long e{0}; e < static_cast<long>(S / sizeof(Tuint)); ++e)
          aln_cache.reductions[e] += overflow_reduction_arr[e];
      }

      for (long v = 0; v < num_vectors; ++v)
      {
        aln_cache.vH_up[v] = simdpp::max(aln_cache.vH_up[v] - overflow_reduction, two_gap_open_pack);
        aln_cache.vF_up[v] = simdpp::max(aln_cache.vF_up[v] - overflow_reduction, two_gap_open_pack);
      }
    }
  }
}


template <typename Tuint>
std::vector<long> inline
get_score_row(AlignmentCache<Tuint> const & aln_cache,
              long const i,
              typename T<Tuint>::vec_pack const & vX)
{
  //AlignmentCache<Tuint> const & aln_cache = *opt.get_alignment_cache();
  //AlignmentResults<Tuint> const & aln_results = *opt.get_alignment_results();
  long const m = aln_cache.query_size;
  long const t = aln_cache.num_vectors;

  assert(t > 0l);
  assert(vX.size() > 0);
  assert(t == static_cast<long>(vX.size()));

  std::vector<Tuint> vec(vX[0].length, 0);
  std::vector<std::vector<Tuint> > mat(t, vec);
  std::vector<long> scores_row;
  scores_row.reserve(m + 1ul);

  for (long v = 0; v < t; ++v)
    simdpp::store_u(&mat[v][0], vX[v]);

  for (long j = 0; j <= m; ++j)
  {
    long const v = j % t;
    long const e = j / t;
    assert(v < static_cast<long>(mat.size()));
    assert(e < static_cast<long>(mat[v].size()));

    long const adjustment = aln_cache.reductions[e] - aln_cache.y_gain * i - aln_cache.x_gain * j;
    scores_row.push_back(static_cast<long>(mat[v][e] + adjustment));
  }

  return scores_row;
}


/*
template <typename Tuint>
void
merge_score_matrices(AlignmentOptions<Tuint> & final_opts, std::vector<AlignmentOptions<Tuint> *> const & opts_vec_ptr)
{
  using Tarr = std::array<long, S / sizeof(Tuint)>;
  using Tpack = typename T<Tuint>::pack;

  assert(opts_vec_ptr.size() >= 1ul);

  final_opts = *opts_vec_ptr[0];
  auto & results = *final_opts.get_alignment_results();
  auto const & cache = *final_opts.get_alignment_cache();

  // Find an array that has the global max of reductions over all input reductions
  Tarr max_reductions;
  max_reductions = opts_vec_ptr[0]->get_alignment_results()->reductions;

  for (long i = 1; i < static_cast<long>(opts_vec_ptr.size()); ++i)
  {
    Tarr const & reductions = opts_vec_ptr[i]->get_alignment_results()->reductions;

    for (long j = 0; j < static_cast<long>(S / sizeof(Tuint)); ++j)
    {
      max_reductions[j] = std::max(reductions[j], max_reductions[j]);
    }
  }


  for (long j = 0; j < static_cast<long>(S / sizeof(Tuint)); ++j)
  {
    std::cerr << j << " " << max_reductions[j] << "\n";
  }


  // Gather how much we need to reduce each vector such they have all reduced the same
  for (long i = 0; i < static_cast<long>(opts_vec_ptr.size()); ++i)
  {
    Tarr const & reductions = opts_vec_ptr[i]->get_alignment_results()->reductions;
    Tarr extra_reductions;
    bool any_extra_reductions = false;

    for (long j = 0; j < static_cast<long>(S / sizeof(Tuint)); ++j)
    {
      extra_reductions[j] = max_reductions[j] - reductions[j];

      if (extra_reductions[j] != 0)
        any_extra_reductions = true;
    }

    if (!any_extra_reductions)
      continue;

    Tpack extra_reductions_vec = simdpp::load(&extra_reductions[0]);
    auto & vF_up = opts_vec_ptr[i]->get_alignment_results()->vF_up;
    auto & vH_up = opts_vec_ptr[i]->get_alignment_results()->vH_up;

    for (long v = 0; v < static_cast<long>(vF_up.size()); ++v)
    {
      vF_up[v] = simdpp::sub_sat(vF_up[v], extra_reductions_vec);
      vH_up[v] = simdpp::sub_sat(vH_up[v], extra_reductions_vec);
    }

    // Set all reductions as "max_reductions"
    opts_vec_ptr[i]->get_alignment_results()->reductions = max_reductions;
  }

  results.reductions = max_reductions;
  results.vF_up.reserve(cache.num_vectors);
  results.vH_up.reserve(cache.num_vectors);

  // Find the best scores at each position
  for (long v = 0; v < cache.num_vectors; ++v)
  {
    assert(opts_vec_ptr[0]);
    assert(opts_vec_ptr[0]->get_alignment_results());
    assert(static_cast<long>(opts_vec_ptr[0]->get_alignment_results()->vF_up.size()) == cache.num_vectors);
    assert(v == static_cast<long>(results.vF_up.size()));
    assert(v == static_cast<long>(results.vH_up.size()));

    results.vF_up.push_back(opts_vec_ptr[0]->get_alignment_results()->vF_up[v]);
    results.vH_up.push_back(opts_vec_ptr[0]->get_alignment_results()->vH_up[v]);

    for (long i = 1; i < static_cast<long>(opts_vec_ptr.size()); ++i)
    {
      results.vF_up[v] = simdpp::max(results.vF_up[v], opts_vec_ptr[i]->get_alignment_results()->vF_up[v]);
      results.vH_up[v] = simdpp::max(results.vH_up[v], opts_vec_ptr[i]->get_alignment_results()->vH_up[v]);
    }
  }
}
*/

#ifndef NDEBUG

template <typename Tuint>
inline void
store_scores(AlignmentOptions<Tuint> & opt,
             AlignmentCache<Tuint> const & aln_cache,
             long const i,
             typename T<Tuint>::vec_pack const & vE)
{
  opt.score_matrix.push_back(get_score_row(aln_cache, i, aln_cache.vH_up));
  opt.vE_scores.push_back(get_score_row(aln_cache, i, vE));
  opt.vF_scores.push_back(get_score_row(aln_cache, i, aln_cache.vF_up));
  assert(opt.score_matrix.size() == opt.score_matrix.size());
  assert(opt.vE_scores.size() == opt.vF_scores.size());
}


#endif // NDEBUG

} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
