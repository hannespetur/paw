#pragma once

#include <cstdint>
#include <type_traits>

#include <simdpp/simd.h>

#include <paw/align/alignment_results.hpp>
#include <paw/align/libsimdpp_utils.hpp>


namespace paw
{

template<typename Tuint>
struct AlignmentOptions
{

public:
  using Tpack = typename T<Tuint>::pack;
  using Tvec_pack = typename T<Tuint>::vec_pack;
  using Tarr_vec_pack = typename T<Tuint>::arr_vec_pack;

  bool left_column_free = false;
  bool right_column_free = false;
  bool continuous_alignment = false; // When set, always continue with the same alignment as long as the query is the same
  std::string query{""};
  long query_size{0}; // Size of query sequence, sometimes also noted as 'm'
  long num_vectors{0}; // Number of SIMD vectors per row
  Tvec_pack vH_up;
  Tvec_pack vF_up;
  Tuint x_gain{0};
  Tuint y_gain{0};
  Tuint top_left_score{0};

  Tuint match_val{0};
  Tuint mismatch_val{0};
  Tuint gap_open_val_x{0};
  Tuint gap_open_val_y{0};
  Tuint gap_open_val{0};
  Tuint max_score_val{0};

  AlignmentResults<Tuint> ar;
  long reduction{0};
  std::array<long, S / sizeof(Tuint)> reductions;
  Tarr_vec_pack W_profile;

private:
  /// User options
  Tuint match = 2;
  Tuint mismatch = 2;
  Tuint gap_open = 5;
  Tuint gap_extend = 1;

  bool backtracking = true;
  bool top_row_free = false;
  bool bottom_row_free = false;
  bool gap_open_free = false;
  bool top_row_gap_open_free = false;
  bool bottom_row_gap_open_free = false;


public:

  AlignmentOptions() = default;


  AlignmentOptions &
  set_match(long val)
  {
    match = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  AlignmentOptions &
  set_mismatch(long val)
  {
    mismatch = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  AlignmentOptions &
  set_gap_open(long val)
  {
    gap_open = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  AlignmentOptions &
  set_gap_extend(long val)
  {
    gap_extend = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  AlignmentOptions &
  set_traceback(bool val)
  {
    backtracking = val;
    return *this;
  }

  inline Tuint get_match() const {return match;}
  inline Tuint get_mismatch() const {return mismatch;}
  inline Tuint get_gap_open() const {return gap_open;}
  inline Tuint get_gap_extend() const {return gap_extend;}
  inline AlignmentResults<Tuint> const & get_results() const {return ar;}


#ifndef NDEBUG
  // Store the score matrix in debug mode
  std::vector<std::vector<long> > score_matrix;
  std::vector<std::vector<long> > vE_scores;
  std::vector<std::vector<long> > vF_scores;
#endif // not NDEBUG

  std::vector<std::vector<long> > merge_solver_vH;
  std::vector<std::vector<long> > merge_solver_vF;

};


namespace SIMDPP_ARCH_NAMESPACE
{


template <typename Tuint>
inline void
init_vH_up(AlignmentOptions<Tuint> & opts)
{
  typename T<Tuint>::vec_uint new_vH0(T<Tuint>::pack::length, 2 * opts.gap_open_val + std::numeric_limits<Tuint>::min());
  assert(opts.vH_up.size() > 0);
  new_vH0[0] = opts.gap_open_val * 3 + std::numeric_limits<Tuint>::min();
  opts.vH_up[0] = simdpp::load_u(&new_vH0[0]);
}


template<typename Tuint, typename Tseq>
void
set_query(AlignmentOptions<Tuint> & opt, Tseq const & seq)
{
  using Tpack = typename T<Tuint>::pack;
  using Tvec_uint = typename T<Tuint>::vec_uint;
  using Tvec_pack = typename T<Tuint>::vec_pack;

  std::string new_query(begin(seq), end(seq));
  Tpack const min_value_pack = simdpp::make_int(std::numeric_limits<Tuint>::min());

  // If it is the same query we can reuse previously calculated numbers
  if (new_query == opt.query)
  {
    if (!opt.continuous_alignment)
    {
      opt.vH_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors),
                            static_cast<Tpack>(simdpp::make_int(2 * opt.gap_open_val + std::numeric_limits<Tuint>::min()))
        );

      init_vH_up(opt);

      opt.vF_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors), min_value_pack);
      opt.reductions.fill(0);
      opt.reduction = 0;

#ifndef NDEBUG
      opt.score_matrix.clear();
      opt.vE_scores.clear();
      opt.vF_scores.clear();
#endif // NDEBUG
    }

    return;
  }

  opt.query = std::move(new_query);
  opt.query_size = opt.query.size();
  opt.num_vectors = (opt.query_size + Tpack::length) / Tpack::length;
  opt.x_gain = opt.get_gap_extend();
  opt.y_gain = std::max(opt.get_gap_extend(),
                                static_cast<Tuint>(opt.get_mismatch() - opt.x_gain));
  opt.gap_open_val_x = opt.get_gap_open() - opt.x_gain;
  opt.gap_open_val_y = opt.get_gap_open() - opt.y_gain;
  opt.gap_open_val = std::max(opt.gap_open_val_x, opt.gap_open_val_y);
  opt.top_left_score = opt.gap_open_val * 3 + std::numeric_limits<Tuint>::min();
  opt.vH_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors),
                        static_cast<Tpack>(simdpp::make_int(2 * opt.gap_open_val + std::numeric_limits<Tuint>::min()))
                        );
  init_vH_up(opt);
  opt.vF_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors), min_value_pack);
  opt.match_val = opt.x_gain + opt.y_gain + opt.get_match();
  opt.mismatch_val = opt.x_gain + opt.y_gain - opt.get_mismatch();
  opt.max_score_val = std::numeric_limits<Tuint>::max() - opt.match_val - opt.gap_open_val;
  opt.reduction = 0;
  opt.reductions.fill(0);


#ifndef NDEBUG
  opt.score_matrix.clear();
  opt.vE_scores.clear();
  opt.vF_scores.clear();
#endif // NDEBUG

  /// Calculate DNA W_profile
  {
    std::array<char, 4> constexpr DNA_BASES = {{'A', 'C', 'G', 'T'}};

    for (std::size_t i = 0; i < DNA_BASES.size(); ++i)
    {
      char const dna_base = DNA_BASES[i];
      auto & W = opt.W_profile[i];
      W.clear(); // Clear previous elements
      W.reserve(opt.num_vectors);

      for (long v = 0; v < opt.num_vectors; ++v)
      {
        Tvec_uint seq(T<Tuint>::pack::length, opt.mismatch_val);

        for (long e = 0, j = v; j < opt.query_size; j += opt.num_vectors, ++e)
        {
          if (dna_base == *(begin(opt.query) + j))
            seq[e] = opt.match_val;
        }

        W.push_back(static_cast<typename T<Tuint>::pack>(simdpp::load_u(&seq[0])));
      }
    }

    assert(static_cast<std::size_t>(opt.num_vectors) == opt.W_profile[0].size());
    assert(static_cast<std::size_t>(opt.num_vectors) == opt.W_profile[1].size());
    assert(static_cast<std::size_t>(opt.num_vectors) == opt.W_profile[2].size());
    assert(static_cast<std::size_t>(opt.num_vectors) == opt.W_profile[3].size());
  } /// Done calculating DNA W_profile
}


template<typename Tuint>
void
reduce_too_high_scores(AlignmentOptions<Tuint> & opt)
{
  using Tpack = typename T<Tuint>::pack;
  using Tmask = typename T<Tuint>::mask;
  using Tarr_uint = typename T<Tuint>::arr_uint;

  long const t = opt.num_vectors;
  Tpack const max_score_pack = simdpp::make_int(opt.max_score_val);

  if (simdpp::reduce_max(opt.vH_up[t - 1]) >= opt.max_score_val)
  {
    /// Reducing value losslessly
    Tarr_uint vF0;
    vF0.fill(0);
    // Store the optimal scores in vector 0
    simdpp::store_u(&vF0[0], opt.vF_up[0]);
    assert(vF0.size() == opt.reductions.size());
    Tarr_uint new_reductions;
    new_reductions.fill(0);
    assert(vF0.size() == new_reductions.size());
    bool any_reductions = false;

    for (long e = 1; e < static_cast<long>(vF0.size()); ++e)
    {
      long new_reduction_val = static_cast<long>(vF0[e]) - static_cast<long>(2 * opt.gap_open_val);

      if (new_reduction_val > 0)
      {
        opt.reductions[e] += new_reduction_val;
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
        opt.vF_up[v] = opt.vF_up[v] - new_reductions_pack;
        opt.vH_up[v] = opt.vH_up[v] - new_reductions_pack;
      }
    }

    Tuint const max_score = simdpp::reduce_max(opt.vH_up[t - 1]);

    if (max_score >= opt.max_score_val)
    {
      Tpack const two_gap_open_pack = simdpp::make_int(2 * opt.gap_open_val);

      /// Reducing values lossily
      Tmask is_about_to_overflow = opt.vH_up[t - 1] >= max_score_pack;

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
          opt.reductions[e] += overflow_reduction_arr[e];
      }

      for (long v = 0; v < t; ++v)
      {
        opt.vH_up[v] = simdpp::max(opt.vH_up[v] - overflow_reduction, two_gap_open_pack);
        opt.vF_up[v] = simdpp::max(opt.vF_up[v] - overflow_reduction, two_gap_open_pack);
      }
    }
  }
}


template<typename Tuint>
std::vector<long> inline
get_score_row(AlignmentOptions<Tuint> const & opt, long const i, typename T<Tuint>::vec_pack const & vX)
{
  using Tvec_uint = typename T<Tuint>::vec_uint;

  long const m = opt.query_size;
  long const t = opt.num_vectors;

  assert(t > 0l);
  assert(t == static_cast<long>(vX.size()));

  Tvec_uint vec(vX[0].length, 0);
  std::vector<Tvec_uint> mat(t, vec);
  std::vector<long> scores_row;
  scores_row.reserve(m + 1ul);

  for (long v = 0; v < t; ++v)
    simdpp::store_u(&mat[v][0], opt.vH_up[v]);

  for (long j = 0; j <= m; ++j)
  {
    long const v = j % t;
    long const e = j / t;
    assert(v < static_cast<long>(mat.size()));
    assert(e < static_cast<long>(mat[v].size()));

    long const adjustment = opt.reductions[e] - opt.top_left_score - opt.y_gain * i - opt.x_gain * j;
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
  std::vector<long> max_vH = get_score_row(final_opts, 0, final_opts.vH_up);
  std::vector<long> max_vF = get_score_row(final_opts, 0, final_opts.vF_up);
  final_opts.merge_solver_vH = std::vector<std::vector<long> >(std::vector<long>(1, 0), max_vH.size());
  final_opts.merge_solver_vF = std::vector<std::vector<long> >(std::vector<long>(1, 0), max_vF.size());

  auto set_max_scores = [&opts_vec_ptr](std::vector<long> & max_scores, std::vector<std::vector<long> > & merge_solver)
  {
    long const in_degree = static_cast<long>(opts_vec_ptr.size()); // InDegree of node that is getting merged
    long const num_scores = static_cast<long>(max_scores.size()); // Number of scores in a row
    merge_solver = std::vector<std::vector<long> >(std::vector<long>(1ul, 0l), max_scores.size());

    for (long i = 1; i < in_degree; ++i)
    {
      AlignmentOptions<Tuint> const & opt = &(opts_vec_ptr[i]);
      std::vector<long> scores = get_score_row(opt, 0, opt.vH_up);
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
    long const max_score_val = final_opts.max_score_val;
    long const max_diff = max_score_val - min_score_val;

    std::array<long, S / sizeof(Tuint)> min_values;
    min_values.fill(std::numeric_limits<long>::max());

    std::array<long, S / sizeof(Tuint)> max_values;
    max_values.fill(std::numeric_limits<long>::min());

    assert(max_vH.size() - 1l == final_opts.query_size);
    assert(max_vF.size() - 1l == final_opts.query_size);

    long t = final_opts.num_vectors;

    for (long j = 0; j <= final_opts.query_size; ++j)
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
  opt.score_matrix.push_back(get_score_row(opt, i, opt.vH_up));
  opt.vE_scores.push_back(get_score_row(opt, i, vE));
  opt.vF_scores.push_back(get_score_row(opt, i, opt.vF_up));
  assert(opt.score_matrix.size() == opt.score_matrix.size());
  assert(opt.vE_scores.size() == opt.vF_scores.size());
}

#endif // NDEBUG

} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
