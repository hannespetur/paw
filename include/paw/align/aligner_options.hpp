#pragma once

namespace paw
{

template <typename Tuint>
struct AlignerOptions
{
  bool default_options; // If set to true, all options are assumed to be the default
  //                       This removes some of the runtime checks

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

  AlignerOptions(bool const _default_options = false)
      : default_options(_default_options)
  {}
};

}