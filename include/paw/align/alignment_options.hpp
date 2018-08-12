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

  bool continuous_alignment = false; // When to true, always continue with the same alignment as long as the query is the same
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
  bool top_row_gap_open_free = gap_open_free;
  bool bottom_row_gap_open_free = gap_open_free;
  bool left_column_gap_open_free = gap_open_free;
  bool right_column_gap_open_free = gap_open_free;

  /// Derived options
  long last_score{0};
  SIMDPP_ARCH_NAMESPACE::Backtrack<Tuint> last_backtrack;


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

  Tuint get_match() const {return match;}
  Tuint get_mismatch() const {return mismatch;}
  Tuint get_gap_open() const {return gap_open;}
  Tuint get_gap_extend() const {return gap_extend;}


#ifndef NDEBUG
  // Store the score matrix in debug mode
  std::vector<std::vector<long> > score_matrix;
#endif // not NDEBUG

};


namespace SIMDPP_ARCH_NAMESPACE
{

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

      opt.vF_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors), min_value_pack);
      opt.reductions.fill(0);

#ifndef NDEBUG
      opt.score_matrix.clear();
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
  opt.vF_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors), min_value_pack);
  opt.match_val = opt.x_gain + opt.y_gain + opt.get_match();
  opt.mismatch_val = opt.x_gain + opt.y_gain - opt.get_mismatch();

  opt.max_score_val = std::numeric_limits<Tuint>::max() - opt.match_val - opt.gap_open_val;
  opt.reductions.fill(0);

#ifndef NDEBUG
  opt.score_matrix.clear();
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

    Tpack const max_scores = opt.vH_up[t - 1];
    Tuint const max_score = simdpp::reduce_max(max_scores);

    if (max_score >= opt.max_score_val)
    {
      Tpack const two_gap_open_pack = simdpp::make_int(2 * opt.gap_open_val);

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


#ifndef NDEBUG

template<typename Tuint>
inline void
store_vH_up_scores(AlignmentOptions<Tuint> & opt,
             long m,
             long const i
  )
{
  using Tvec_uint = typename T<Tuint>::vec_uint;

  long const t = opt.vH_up.size();
  assert(t > 0);
  Tvec_uint vec(opt.vH_up[0].length, 0);
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

    scores_row.push_back(static_cast<long>(mat[v][e] - opt.top_left_score - opt.y_gain * i - opt.x_gain * j + opt.reductions[e]));
  }

  //for (auto const val : scores_row)
  //  std::cout << std::setw(5) << val;
  //std::cout << std::endl;
  opt.score_matrix.push_back(std::move(scores_row));
}

#endif // NDEBUG

} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
