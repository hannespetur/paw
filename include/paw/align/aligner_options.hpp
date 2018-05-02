#pragma once

#include <cstdint>

#include <paw/align/libsimdpp_utils.hpp>


namespace paw
{

struct AlignerOptions
{
  #if defined(PAW_USE_UINT8)
  using uint = uint8_t;
  #else
  using uint = uint16_t;
  #endif

  bool default_options = false; // If set to true, all options are assumed to be the default
                                // This can remove some of the runtime checks

  uint match = 2;
  uint mismatch = 2;
  uint gap_open = 5;
  uint gap_extend = 1;

  bool backtracking = true;
  bool top_row_free = false;
  bool bottom_row_free = false;
  bool gap_open_free = false;
  bool top_row_gap_open_free = gap_open_free;
  bool bottom_row_gap_open_free = gap_open_free;
  bool left_column_gap_open_free = gap_open_free;
  bool right_column_gap_open_free = gap_open_free;

  explicit AlignerOptions(bool const _default_options = false)
    : default_options(_default_options)
  {}


  void
  set_match(int val)
  {
    match = val >= 0 ? val : -val;
  }


  void
  set_mismatch(int val)
  {
    mismatch = val >= 0 ? val : -val;
  }


  void
  set_gap_open(int val)
  {
    gap_open = val >= 0 ? val : -val;
  }


  void
  set_gap_extend(int val)
  {
    gap_extend = val >= 0 ? val : -val;
  }


};

} // namespace paw
