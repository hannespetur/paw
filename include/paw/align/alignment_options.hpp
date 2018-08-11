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
  std::string query{""};
  long query_size{0}; // Size of query sequence

private:
  /// User options
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

  /// Derived options
  long last_score{0};
  SIMDPP_ARCH_NAMESPACE::Backtrack<Tuint> last_backtrack;
  typename T<Tuint>::vec_pack last_vH;
  typename T<Tuint>::vec_pack last_vF;


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


  template<typename Tseq>
  void
  set_query(Tseq const & seq)
  {
    {
      std::string new_query(begin(seq), end(seq));

      // If it is the same query we can reuse previously calculated numbers
      if (new_query == query)
        return;

      query = std::move(new_query);
    }

    query_size = query.size();
  }


  AlignmentOptions &
  set_traceback(bool val)
  {
    backtracking = val;
    return *this;
  }

  uint get_match() const {return match;}
  uint get_mismatch() const {return mismatch;}
  uint get_gap_open() const {return gap_open;}
  uint get_gap_extend() const {return gap_extend;}

};

} // namespace paw
