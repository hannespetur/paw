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

  std::string query{""};
  long query_size{0}; // Size of query sequence, sometimes also noted as 'm'
  long num_vectors{0}; // Number of SIMD vectors per row
  Tvec_pack vH_up;
  Tvec_pack vF_up;
  Tuint x_gain{0};
  Tuint y_gain{0};

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
  set_match(int val)
  {
    match = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  AlignmentOptions &
  set_mismatch(int val)
  {
    mismatch = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  AlignmentOptions &
  set_gap_open(int val)
  {
    gap_open = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
    return *this;
  }


  AlignmentOptions &
  set_gap_extend(int val)
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


//#ifndef NDEBUG
  // Store the score matrix in debug mode
  std::vector<std::vector<long> > score_matrix;
//#endif // not NDEBUG

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

  // If it is the same query we can reuse previously calculated numbers
  if (new_query == opt.query)
    return;

  Tpack const min_value_pack = simdpp::make_int(std::numeric_limits<Tuint>::min());
  opt.query = std::move(new_query);
  opt.query_size = opt.query.size();
  opt.num_vectors = (opt.query_size + Tpack::length) / Tpack::length;
  opt.x_gain = opt.get_gap_extend();
  opt.y_gain = std::max(opt.get_gap_extend(),
                                static_cast<Tuint>(opt.get_mismatch() - opt.x_gain));
  opt.gap_open_val_x = opt.get_gap_open() - opt.x_gain;
  opt.gap_open_val_y = opt.get_gap_open() - opt.y_gain;
  opt.gap_open_val = std::max(opt.gap_open_val_x, opt.gap_open_val_y);
  opt.vH_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors),
                        static_cast<Tpack>(simdpp::make_int(2 * opt.gap_open_val + std::numeric_limits<Tuint>::min()))
                        );
  opt.vF_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors), min_value_pack);
  opt.match_val = opt.x_gain + opt.y_gain + opt.get_match();
  opt.mismatch_val = opt.x_gain + opt.y_gain - opt.get_mismatch();

  opt.max_score_val = std::numeric_limits<Tuint>::max() - opt.match_val - opt.gap_open_val;
  long const top_left_score = opt.gap_open_val * 3 + std::numeric_limits<Tuint>::min();
  opt.reductions.fill(-top_left_score);

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


#ifndef NDEBUG

template<typename Tuint>
inline void
store_scores(AlignmentOptions<Tuint> & opt,
             long m,
             typename T<Tuint>::vec_pack const & vX,
             long const i
  )
{
  using Tvec_uint = typename T<Tuint>::vec_uint;

  long const t = vX.size();
  assert(t > 0);
  Tvec_uint vec(vX[0].length, 0);
  std::vector<Tvec_uint> mat(t, vec);
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

    scores_row.push_back(static_cast<long>(mat[v][e] - opt.y_gain * i - opt.x_gain * j + opt.reductions[e]));
  }

  //for (auto const val : scores_row)
  //  std::cout << std::setw(5) << val;
  //std::cout << std::endl;
  opt.score_matrix.push_back(std::move(scores_row));
}

#endif // NDEBUG

} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
