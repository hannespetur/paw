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
  using uint = typename std::make_unsigned<Tuint>::type;

  bool default_options = false; // If set to true, all options are assumed to be the default
                                // This can remove some of the runtime checks

  bool backtracking = true;
  bool top_row_free = false;
  bool bottom_row_free = false;
  bool gap_open_free = false;
  bool top_row_gap_open_free = gap_open_free;
  bool bottom_row_gap_open_free = gap_open_free;
  bool left_column_gap_open_free = gap_open_free;
  bool right_column_gap_open_free = gap_open_free;


private:
  uint match = 2;
  uint mismatch = 2;
  uint gap_open = 5;
  uint gap_extend = 1;

  long last_score{0};
  SIMDPP_ARCH_NAMESPACE::Backtrack<Tuint> last_backtrack;
  typename T<Tuint>::vec_pack last_vH;
  typename T<Tuint>::vec_pack last_vF;


public:

  explicit AlignmentOptions(bool const _default_options = false)
    : default_options(_default_options)
  {}


  void
  set_match(int val)
  {
    match = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
  }


  void
  set_mismatch(int val)
  {
    mismatch = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
  }


  void
  set_gap_open(int val)
  {
    gap_open = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
  }


  void
  set_gap_extend(int val)
  {
    gap_extend = val >= 0 ? static_cast<Tuint>(val) : static_cast<Tuint>(-val);
  }

  uint get_match() const {return match;}
  uint get_mismatch() const {return mismatch;}
  uint get_gap_open() const {return gap_open;}
  uint get_gap_extend() const {return gap_extend;}

};

} // namespace paw
