#pragma once

#include <array> // std::array
#include <cassert>
#include <cstdint> // uint8_t
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>


#include <paw/align/aligner_options.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>

#include <simdpp/simd.h>


namespace paw
{

//namespace simd = simdpp;


struct Row2
{
  using Tint = uint16_t;
  //using Tpack = simdpp::int16<16 / sizeof(int16_t), void>;
  //using Tmask = simdpp::mask_int16<16 / sizeof(int16_t), void>;
  using Tpack = simdpp::int16<16 / sizeof(Tint), void>;
  using Tlogical_pack = simdpp::mask_int16<16 / sizeof(Tint), void>;
  using Tvec = std::vector<Tpack, simdpp::aligned_allocator<Tpack, sizeof(Tpack)> >;

  std::size_t static const vector_size = Tpack::length;

  std::size_t const n_elements = 0;
  Tvec vectors;

  /* CONSTRUCTORS */
  Row2() = default;

  Row2(std::size_t const _n_elements)
    : n_elements(_n_elements)
  {
    Tpack my_vector = simdpp::make_zero();
    vectors.resize((n_elements + vector_size - 1) / vector_size, my_vector);
  }

  Row2(std::size_t const _n_elements, Tint const val)
    : n_elements(_n_elements)
  {
    Tpack my_vector = simdpp::make_int(val);
    vectors.resize((n_elements + vector_size - 1) / vector_size, my_vector);
  }
};



template <typename Tit>
class Align
{
public:
  using Trow = Row2; // A row of vectors that can be run in parallel
  using Tarr = std::array<Trow, 4>;
  using Tpack = typename Trow::Tpack;
  using Tlogical_pack = typename Trow::Tlogical_pack;
  using Tint = typename Trow::Tint;
  using Tuint = typename Trow::Tint;

  // p is the length (or cardinality) of the SIMD vectors
  std::size_t static const p = Tpack::length;
  void calculate_DNA_W_profile();

  inline Tarr const & get_W_profile() const
  {
    return this->W_profile;
  }


private:
  AlignerOptions<Tuint> const opt;
  std::vector<Event> free_snp_edits;
  std::vector<Event> free_del_edits;
  std::vector<Event> free_ins_edits;
  Tit d_begin;
  Tit d_end;
  Tit q_begin;
  Tit q_end;

  std::size_t const m = 0;
  std::size_t alignment_end = m;
  std::size_t n = 0;
  std::size_t const t = 0;

  Tuint const x_gain;
  Tpack const x_gain_pack = simdpp::make_zero();
  Tuint const y_gain;
  Tpack const y_gain_pack = simdpp::make_zero();

  Tuint const gap_open_val_x = 0;
  Tpack const gap_open_pack_x = simdpp::make_zero();
  Tuint const gap_open_val_y = 0;
  Tpack const gap_open_pack_y = simdpp::make_zero();
  Tuint const gap_open_val = 0;
  Tpack const gap_open_pack = simdpp::make_zero();

  Tuint const max_score_val = 0;
  Tpack const max_score_pack = simdpp::make_zero();

  int64_t total_reductions = 0;
  Trow vH_up; // Previous H row
  Trow vH;    // Current H row
  Trow vE;    // Current E row
  Trow vF_up; // Previous F row
  Trow vF;    // Current F row
  Tarr W_profile;
  Backtrack<Tuint> mB; //(n /*n_row*/, t /*n_vectors in each score row*/);
  Tuint top_left_score = 0;

  Trow const & get_W(const char b);

  Tpack max_greater(Tpack & v1, Tpack const & v2);
  // Tpack max_greater_or_equal(Tpack & v1, Tpack const & v2);
  void calculate_scores();

  // Same as 'calculate_scores', but assumes default options (less runtime checks)
  // void calculate_scores_default();
  void check_gap_extend_deletions();
  void check_gap_extend_deletions_with_backtracking(std::size_t const i);
  void init_score_vectors();

public:
  Align(Tit _d_begin,
        Tit _d_end,
        AlignerOptions<Tuint> const & _opt = AlignerOptions<Tuint>(true /*default options*/),
        std::set<Event> const & /*free_edits*/ = std::set<Event>()
        )
    : opt(_opt)
    , d_begin(_d_begin)
    , d_end(_d_end)
    , m(std::distance(d_begin, d_end))
    , n(0)
    , t((m + p) / p)
    , x_gain(opt.gap_extend)
    , x_gain_pack(simdpp::make_int(x_gain))
    , y_gain(std::max(opt.gap_extend, static_cast<Tuint>(opt.mismatch - x_gain)))
    , y_gain_pack(simdpp::make_int(y_gain))
    , gap_open_val_x(opt.gap_open - x_gain)
    , gap_open_pack_x(simdpp::make_int(gap_open_val_x))
    , gap_open_val_y(opt.gap_open - y_gain)
    , gap_open_pack_y(simdpp::make_int(gap_open_val_y))
    , gap_open_val(std::max(gap_open_val_x, gap_open_val_y))
    , gap_open_pack(simdpp::make_int(gap_open_val))
    , max_score_val(std::numeric_limits<Tuint>::max() - gap_open_val)
    , max_score_pack(simdpp::make_int(max_score_val))
    , vH_up(m + 1)
    , vH(m + 1)
    , vE(m + 1)
    , vF_up(m + 1)
    , vF(m + 1)
    , W_profile
    {
      {
        Trow(m + 1),
        Trow(m + 1),
        Trow(m + 1),
        Trow(m + 1)
      }
    }
    , mB()
  {
    /*
    for (auto const & e : free_edits)
    {
      if (e.ref.size() == 0)
      {
        free_ins_edits.push_back(e);
      }
      else if (e.alt.size() == 0)
      {
        free_del_edits.push_back(e);
      }
      else   //if (e.ref.size() == 1 && e.alt.size() == 1)
      {
        free_snp_edits.push_back(e);
      }
    }

    auto get_edit_str = [](Event const & e) -> std::string
                        {
                          std::stringstream ss;
                          ss << e.pos << " "
                             << (e.ref.size() > 0 ? e.ref : "-") << " "
                             << (e.alt.size() > 0 ? e.alt : "-");
                          return ss.str();
                        };


    if (free_snp_edits.size() > 0)
    {
      std::cout << "Free SNP edits are:\n";

      for (auto const & e : free_snp_edits)
        std::cout << get_edit_str(e) << "\n";
    }

    if (free_ins_edits.size() > 0)
    {
      std::cout << "Free INS edits are:\n";

      for (auto const & e : free_ins_edits)
        std::cout << get_edit_str(e) << "\n";
    }

    if (free_del_edits.size() > 0)
    {
      std::cout << "Free DEL edits are:\n";

      for (auto const & e : free_del_edits)
        std::cout << get_edit_str(e) << "\n";
    }
    */
  }

  int64_t align(Tit q_begin, Tit q_end);   // Align query
  std::pair<std::string, std::string> get_aligned_strings();
  std::set<Event> get_edit_script(std::pair<std::string, std::string> const & s);
  std::vector<Cigar> get_cigar(std::pair<std::string, std::string> const & s);
};


} // namespace paw

namespace
{

template <typename Tpack>
inline void
print_score_vector_standard(Tpack const & vX)
{
  std::size_t const t = vX.vectors.size();

  for (std::size_t j = 0; j < vX.n_elements; ++j)
  {
    std::size_t const v = j % t;
    std::size_t const e = j / t;
    std::cout << std::setw(4) << static_cast<uint32_t>(vX.vectors.at(v)[e]) << " ";
  }

  std::cout << "\n";
}


template <typename Tpack>
inline void
print_score_vector_vectorized(Tpack const & v)
{
  for (auto const & vec : v.vectors)
    std::cout << vec << " ";

  std::cout << "\n";
}


template <typename Tpack>
inline void
print_score_vectors(Tpack const & vH,
                    Tpack const & vH_up,
                    Tpack const & vE,
                    Tpack const & vF,
                    Tpack const & vF_up,
                    Tpack const & vW
                    )
{
  std::cout << "Standard H_up  : "; print_score_vector_standard(vH_up);
  std::cout << "Standard H     : "; print_score_vector_standard(vH);
  std::cout << "Standard E     : "; print_score_vector_standard(vE);
  std::cout << "Standard F_up  : "; print_score_vector_standard(vF_up);
  std::cout << "Standard F     : "; print_score_vector_standard(vF);
  std::cout << "Standard W     : "; print_score_vector_standard(vW);
  std::cout << "\n";
  std::cout << "Vectorized H_up: "; print_score_vector_vectorized(vH_up);
  std::cout << "Vectorized H   : "; print_score_vector_vectorized(vH);
  std::cout << "Vectorized E   : "; print_score_vector_vectorized(vE);
  std::cout << "=====\n";
}


} // anon namespace


namespace paw
{

template <typename Tit>
Row2 const &
Align<Tit>::get_W(const char b)
{
  switch (b)
  {
  case 'A':
    return W_profile[0];

  case 'C':
    return W_profile[1];

  case 'G':
    return W_profile[2];

  case 'T':
    return W_profile[3];

  default:
    std::cerr << "Warning: " << b << " is not A, C, G, or T.\n";
    assert(false);
    return W_profile[3];
  }
}


template <typename Tit>
typename Align<Tit>::Tpack inline
Align<Tit>::max_greater(Tpack & v1, Tpack const & v2)
{
  Tpack is_greater_pack = v2 > v1;
  Tlogical_pack is_greater = simdpp::to_mask(is_greater_pack);
  v1 = simdpp::blend(v2, v1, is_greater);
  //v1 = boost::simd::if_else(is_greater, v2, v1);
  return simdpp::blend(static_cast<Tpack>(simdpp::make_int(1)),
                       static_cast<Tpack>(simdpp::make_zero()),
                       is_greater
                       );
  //return boost::simd::if_one_else_zero(is_greater);
}


template <typename Tit>
void inline
Align<Tit>::calculate_DNA_W_profile()
{
  std::array<char, 4> constexpr DNA_BASES = {{'A', 'C', 'G', 'T'}};
  Tuint const match = x_gain + y_gain + opt.match;
  Tuint const mismatch = x_gain + y_gain - opt.mismatch;
  long const t = W_profile[0].vectors.size();
  assert(t > 0);

  for (std::size_t i = 0; i < DNA_BASES.size(); ++i)
  {
    char const dna_base = DNA_BASES[i];
    auto & W = W_profile[i];

    for (long v = 0; v < t; ++v)
    {
      std::vector<Tint, simdpp::aligned_allocator<Tint, sizeof(Tint)> > seq;
      seq.reserve(p);

      for (long j = v; j < static_cast<long>(m); j += t)
      {
        if (dna_base == *(d_begin + j))
          seq.push_back(match);
        else
          seq.push_back(mismatch);
      }

      seq.resize(p, mismatch); // Make sure the vector is fully extended
      W.vectors[v] = simdpp::load(&seq[0]);
    }

    // Update the W_profile if there are any free mismatches
    /*
    for (auto const & snp_e : free_snp_edits)
    {
      assert(snp_e.ref.size() == 1);
      assert(snp_e.alt.size() == 1);
      assert(snp_e.pos < m);
      assert(snp_e.ref[0] == *std::next(d_begin, snp_e.pos));

      if (a != snp_e.alt[0])
        continue;

      std::size_t const v = snp_e.pos % t;
      std::size_t const e = snp_e.pos / t;

      W.vectors[v][e] = x_gain + y_gain + match;
    }
    */
  }
}


template <typename Tit>
inline void
Align<Tit>::init_score_vectors()
{
  // For less verbosity of code
  auto const x = gap_open_val;

  if (opt.top_row_gap_open_free)
  {
    for (auto & v : vH_up.vectors)
      v = simdpp::make_int(x * 2);

    for (auto & v : vH.vectors)
      v = simdpp::make_int(x * 2);
  }
  else
  {
    for (auto & v : vH_up.vectors)
      v = simdpp::make_int(x * 2, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x);

    for (auto & v : vH.vectors)
      v = simdpp::make_int(x);
  }

  for (auto & v : vF_up.vectors)
    v = simdpp::make_zero();

  // TODO: Implement
  //if (opt.top_row_free)
  //{
  //  for (std::size_t v = 0; v < t; ++v)
  //  {
  //    for (std::size_t e = 0, j = v; j < m; ++e, j += t)
  //    {
  //      Tuint const gain = j * x_gain;
  //      vH_up.vectors[v][e] = 2 * gap_open_val + gain;
  //      vH.vectors[v][e] = 2 * gap_open_val + gain;
  //      vF_up.vectors[v][e] = gap_open_val + gain;
  //    }
  //  }
  //}

  top_left_score = x * 2;
}


template <typename Tit>
int64_t inline
Align<Tit>::align(Tit _q_begin, Tit _q_end)
{
  q_begin = _q_begin;
  q_end = _q_end;
  n = std::distance(q_begin, q_end);
  this->init_score_vectors();
  total_reductions = 0;

  if (opt.backtracking)
  {
    if (n > mB.matrix.size())
    {
      mB = Backtrack<Tuint>(n, t);
    }
    else
    {
      // Make sure the first n rows have only zeros
      for (std::size_t i = 0; i < n; ++i)
      {
        for (auto & v : mB.matrix[i])
          v = simdpp::make_zero();
      }
    }
  }

  calculate_scores();

  return 0;
}


template <typename Tit>
void inline
Align<Tit>::calculate_scores()
{
  bool const backtracking = opt.backtracking;
  bool const left_column_gap_open_free = opt.left_column_gap_open_free;
  bool const right_column_gap_open_free = opt.right_column_gap_open_free;
  std::size_t const right_v = m % t; // Vector that contains the rightmost element
  std::size_t const right_e = m / t; // The right-most element (in vector 'right_v')
  Tuint const max_gain_per_row = std::max(static_cast<Tuint>(0),
                                          static_cast<Tuint>(opt.match + x_gain + y_gain)
                                          );

  //Tuint current_max_score = vH_up.vectors[0][0]; // Keep track of current max score
  Tuint current_max_score = gap_open_val * 2;

  /// Start of outer loop
  for (std::size_t i = 0; i < n; ++i)
  {
    // vW_i,j has the scores for each substitution between bases q[i] and d[j]
    Trow const & vW = this->get_W(*std::next(q_begin, i));

    // Clear vE vectors
    for (auto & vE_vec : vE.vectors)
      vE_vec = simdpp::make_zero();

    /// Calculate vector 0
    vH.vectors[0] = simdpp::move8_r<1>(vH_up.vectors[t - 1] + vW.vectors[t - 1]);

    /*
    // This is the only way to reach the left-most coordinate
    vH.vectors[0][0] = std::max(static_cast<Tuint>(vF_up.vectors[0][0]),
                                static_cast<Tuint>(vH_up.vectors[0][0] - gap_open_val)
                                );

    // Check if any insertion have highest values
    if (backtracking)
    {
      mB.matrix[i][0][0] = mB.INS_BT;
      vF.vectors[0] = vH_up.vectors[0] - gap_open_pack;

      if (left_column_gap_open_free)
        vF.vectors[0][0] = vH_up.vectors[0][0];

      if (right_column_gap_open_free && right_v == 0)
        vF.vectors[0][right_e] = vH_up.vectors[0][right_e];

      mB.set_ins_extend(i, 0, max_greater(vF.vectors[0], vF_up.vectors[0]));
      mB.set_ins(i, 0, max_greater(vH.vectors[0], vF.vectors[0]));
    }
    else
    {
      vF.vectors[0] = vH_up.vectors[0] - gap_open_pack;

      if (left_column_gap_open_free)
        vF.vectors[0][0] = vH_up.vectors[0][0];

      if (right_column_gap_open_free && right_v == 0)
        vF.vectors[0][right_e] = vH_up.vectors[0][right_e];

      vF.vectors[0] = boost::simd::max(vF.vectors[0], vF_up.vectors[0]);
      vH.vectors[0] = boost::simd::max(vH.vectors[0], vF.vectors[0]);
    }
    /// Done calculating vector 0

    /// Calculate vectors 1,...,v-1
    for (std::size_t v = 1; v < t; ++v)
    {
      // Check for substitutions and if it has a higher score than the insertion
      vH.vectors[v] = vH_up.vectors[v - 1] + vW.vectors[v - 1];
      //vF.vectors[v] = boost::simd::max(vH_up.vectors[v] - gap_open_pack, vF_up.vectors[v]);

      if (backtracking)
      {
        vF.vectors[v] = vH_up.vectors[v] - gap_open_pack;

        // Slight performance gain if this if statement is moved out of the loop
        if (right_column_gap_open_free && v == right_v)
          vF.vectors[right_v][right_e] = vH_up.vectors[right_v][right_e];

        mB.set_ins_extend(i, v, max_greater(vF.vectors[v], vF_up.vectors[v]));
        mB.set_ins(i, v, max_greater(vH.vectors[v], vF.vectors[v]));
      }
      else
      {
        vF.vectors[v] = boost::simd::max(vH_up.vectors[v] - gap_open_pack, vF_up.vectors[v]);

        // Slight performance gain if this if statement is moved out of the loop
        if (right_column_gap_open_free && v == right_v)
          vF.vectors[right_v][right_e] = vH_up.vectors[right_v][right_e];

        vH.vectors[v] = boost::simd::max(vH.vectors[v], vF.vectors[v]);
      }

      // Deletions pass 1: Gap opens
      vE.vectors[v] = vH.vectors[v - 1] - gap_open_pack;
    } /// Done calculating vectors 1,...,v-1

    // Calculate vE.vector[0]
    vE.vectors[0] =
      boost::simd::shuffle<boost::simd::pattern<shift_elements_right> >(
        vH.vectors[t - 1] - gap_open_pack_x);

    //for (std::size_t e = 1; e < p; ++e)
    //  vE.vectors[0][e] = vH.vectors[t - 1][e - 1] - gap_open_val_x;

    // Deletions pass 2: Gap extends
    if (backtracking)
    {
      check_gap_extend_deletions_with_backtracking(i);

      // Check if any vE has higher scores than vH
      for (std::size_t v = 0; v < t; ++v)
        mB.set_del(i, v, max_greater(vH.vectors[v], vE.vectors[v]));
    }
    else
    {
      check_gap_extend_deletions();

      // Check if any vE has higher scores than vH
      for (std::size_t v = 0; v < t; ++v)
        vH.vectors[v] = boost::simd::max(vH.vectors[v], vE.vectors[v]);
    }

    //std::cout << "q[" << i << "] = " << q[i] << "\n";
    //print_score_vectors(vH, vH_up, vE, vF, vF_up, vW); // Useful when debugging
    //std::cout << "max score = " << static_cast<uint16_t>(current_max_score) << "\n";

    if (current_max_score + max_gain_per_row < max_score_val)
    {
      // Assume the worst-case and only find the max score when it would be possible to overflow
      current_max_score += max_gain_per_row;
    }
    else
    {
      // Find the correct max score (accurate but computationally expensive)
      current_max_score = 0;

      for (auto const & vec : vH.vectors)
      {
        current_max_score = std::max(current_max_score, boost::simd::maximum(vec));

        if (current_max_score >= max_score_val)
        {
          current_max_score -= gap_open_val;
          total_reductions += gap_open_val;

          for (auto & vH_vec : vH.vectors)
            vH_vec = boost::simd::max(vH_vec - gap_open_pack, gap_open_pack);

          for (auto & vF_vec : vF.vectors)
            vF_vec = boost::simd::max(vF_vec - gap_open_pack, gap_open_pack);

          break;
        }
      }
    }
    */

    std::swap(vF.vectors, vF_up.vectors);
    std::swap(vH.vectors, vH_up.vectors);
  } /// End of outer loop

}


} // namespace paw
