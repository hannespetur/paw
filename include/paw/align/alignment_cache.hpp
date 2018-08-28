#pragma once

#include <string>

#include <paw/align/libsimdpp_utils.hpp>


namespace paw
{

template<typename Tuint>
struct AlignmentCache
{
  using Tarr_vec_pack = typename T<Tuint>::arr_vec_pack;

  std::string query{""};
  Tuint x_gain{0};
  Tuint y_gain{0};
  long query_size{0}; // Size of query sequence, sometimes also noted as 'm'
  long num_vectors{0}; // Number of SIMD vectors per row
  Tuint match_val{0};
  Tuint mismatch_val{0};
  Tuint gap_open_val{0};
  Tuint max_score_val{0};
  Tarr_vec_pack W_profile;


  inline void
  set_query(std::string && new_query)
  {
    query = std::move(new_query);
    query_size = query.size();
    num_vectors = (query_size + T<Tuint>::pack::length) / T<Tuint>::pack::length;
  }

  inline void
  set_options(Tuint match, Tuint mismatch, Tuint gap_open, Tuint gap_extend)
  {
    x_gain = gap_extend;
    y_gain = std::max(gap_extend, static_cast<Tuint>(mismatch - x_gain));
    gap_open_val = std::max(gap_open - x_gain, gap_open - y_gain);
    match_val = x_gain + y_gain + match;
    mismatch_val = x_gain + y_gain - mismatch;
    max_score_val = std::numeric_limits<Tuint>::max() - match_val - gap_open_val;

    /// Calculate DNA W_profile
    {
      std::array<char, 4> constexpr DNA_BASES = {{'A', 'C', 'G', 'T'}};

      for (std::size_t i = 0; i < DNA_BASES.size(); ++i)
      {
        char const dna_base = DNA_BASES[i];
        auto & W = W_profile[i];
        W.clear(); // Clear previous elements
        W.reserve(num_vectors);

        for (long v = 0; v < num_vectors; ++v)
        {
          typename T<Tuint>::vec_uint seq(T<Tuint>::pack::length, mismatch_val);

          for (long e = 0, j = v; j < query_size; j += num_vectors, ++e)
          {
            if (dna_base == *(begin(query) + j))
              seq[e] = match_val;
          }

          W.push_back(static_cast<typename T<Tuint>::pack>(simdpp::load_u(&seq[0])));
        }
      }

      assert(static_cast<std::size_t>(num_vectors) == W_profile[0].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[1].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[2].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[3].size());
    } /// Done calculating DNA W_profile
  }

};

} // namespace paw