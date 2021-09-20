#pragma once

#include <string>

// cppcheck-suppress preprocessorErrorDirective
#if defined __has_include
// cppcheck-suppress preprocessorErrorDirective
#  if __has_include (<string_view>)
#    include <string_view>
#  endif
#endif

#include <paw/align/event.hpp>
#include <paw/align/libsimdpp_utils.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>

namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

template <typename Tuint>
struct AlignmentCache
{
  using Tarr_vec_pack = typename T<Tuint>::arr_vec_pack;
  using Tarr_uint = typename T<Tuint>::arr_uint;
  using Tvec_pack = typename T<Tuint>::vec_pack;

  std::string query;

  Tuint x_gain {0};
  Tuint y_gain {0};
  long query_size {0}; // Size of query sequence, sometimes also noted as 'm'
  long num_vectors {0}; // Number of SIMD vectors per row
  Tuint match_val {0};
  Tuint mismatch_val {0};
  Tuint gap_open_val {0};
  Tuint max_score_val {0};
  Tvec_pack vH_up{};
  Tvec_pack vF_up{};
  Backtrack<Tuint> mB{};
  Tarr_vec_pack W_profile;
  std::array<long, S / sizeof(Tuint)> reductions;

  AlignmentCache()
  {
    //W_profile.fill(0);
    reductions.fill(0);
  }


  inline void
  set_query(std::string && new_query)
  {
    query = std::forward<std::string>(new_query);
    query_size = query.size();
    num_vectors = (query_size + T<Tuint>::pack::length) /
                  T<Tuint>::pack::length;
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

      for (long i = 0; i < 4; ++i)
      {
        char const dna_base = DNA_BASES[i];
        auto & W = W_profile[i];
        W.clear(); // Clear previous elements
        W.reserve(num_vectors);

        {
          for (long v = 0; v < num_vectors; ++v)
          {
            std::vector<Tuint> seq(T<Tuint>::pack::length, mismatch_val);

            for (long e = 0, j = v; j < query_size; j += num_vectors, ++e)
            {
              assert(j < static_cast<long>(query.size()));
              char const query_dna_base = *(begin(query) + j);

              if (dna_base == query_dna_base || query_dna_base == 'N')
                seq[e] = match_val;
            }

            W.push_back(static_cast<typename T<Tuint>::pack>(simdpp::load_u(&seq[0])));
          }
        }
      }

      // All is a match with N
      {
        auto & W = W_profile[4];
        W.clear(); // Clear previous elements
        W.reserve(num_vectors);

        for (long v{0}; v < num_vectors; ++v)
        {
          std::vector<Tuint> seq(T<Tuint>::pack::length, match_val);
          W.push_back(static_cast<typename T<Tuint>::pack>(simdpp::load_u(&seq[0])));
        }
      }

      // All is a mismatch with X
      {
        auto & W = W_profile[5];
        W.clear(); // Clear previous elements
        W.reserve(num_vectors);

        for (long v{0}; v < num_vectors; ++v)
        {
          std::vector<Tuint> seq(T<Tuint>::pack::length, mismatch_val);
          W.push_back(static_cast<typename T<Tuint>::pack>(simdpp::load_u(&seq[0])));
        }
      }

      assert(static_cast<std::size_t>(num_vectors) == W_profile[0].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[1].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[2].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[3].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[4].size());
      assert(static_cast<std::size_t>(num_vectors) == W_profile[5].size());
    } /// Done calculating DNA W_profile
  }


  inline void
  reduce_every_element(long val)
  {
    std::for_each(reductions.begin(), reductions.end(), [val](long & element){element += val;});
  }


  inline void
  set_free_snps(std::set<Event2> const & events)
  {
    for (auto const & e : events)
      set_free_snp(e.pos, e.alt[0]);
  }


  inline void
  set_free_snp(long pos, char alt)
  {
    assert(pos < static_cast<long>(query.size()));
    //std::cerr << "Adding free SNP pos,alt = " << pos << ',' << alt << '\n';
    auto & W = W_profile[magic_function(alt)];
    assert(pos < static_cast<long>(query.size()));

    long const v = pos % num_vectors;
    long const e = pos / num_vectors;

    Tarr_uint vW;
    vW.fill(std::numeric_limits<Tuint>::min());
    simdpp::store_u(&vW[0], W[v]);
    //std::cout << "vW[e] changed from " << static_cast<int>(vW[e]) << " to " << static_cast<int>(match_val) << "\n";
    vW[e] = match_val + 1;
    W[v] = simdpp::load_u(&vW[0]);
  }


};

} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
