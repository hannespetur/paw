#pragma once

#include <cstdint>
#include <type_traits>


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

  Tpack gap_open_pack_x = simdpp::make_int(0);
  Tpack gap_open_pack_y = simdpp::make_int(0);
  Tpack const min_value_pack = simdpp::make_int(std::numeric_limits<Tuint>::min());

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

};


namespace SIMDPP_ARCH_NAMESPACE
{

template<typename Tuint, typename Tseq>
void
set_query(AlignmentOptions<Tuint> & opt, Tseq const & seq)
{
  using Tpack = typename T<Tuint>::pack;
  using Tvec_pack = typename T<Tuint>::vec_pack;

  std::string new_query(begin(seq), end(seq));

  // If it is the same query we can reuse previously calculated numbers
  if (new_query == opt.query)
    return;

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
  opt.vF_up = Tvec_pack(static_cast<std::size_t>(opt.num_vectors), opt.min_value_pack);
  opt.gap_open_pack_x = simdpp::make_int(opt.gap_open_val_x);
  opt.gap_open_pack_y = simdpp::make_int(opt.gap_open_val_y);
  opt.match_val = opt.x_gain + opt.y_gain + opt.get_match();
  opt.mismatch_val = opt.x_gain + opt.y_gain - opt.get_mismatch();
}

} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw
