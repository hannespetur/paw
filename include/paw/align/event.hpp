#pragma once

#include <algorithm>
#include <cassert>
#include <set>
#include <string> // std::string
#include <vector>

namespace paw
{

// An event is an single point mutation with two alleles, reference and alternative.
class Event2
{
public:
  long pos{0}; // Position of the event, compared to the database.
  long pos_q{0}; // Position of the event, compared to the query.
  std::string ref; // Reference allele.
  std::string alt; // Alternative allele.

  Event2(long _pos, long _pos_q, std::string const & _ref, std::string const & _alt)
    : pos(_pos)
    , pos_q(_pos_q)
    , ref(_ref)
    , alt(_alt)
  {}

  inline bool
  is_snp() const
  {
    return ref.size() == 1 && alt.size() == 1;
  }


  inline bool
  is_deletion() const
  {
    return alt.size() == 0;
  }


  inline bool
  is_insertion() const
  {
    return ref.size() == 0;
  }

  inline bool
  is_indel() const
  {
    return is_deletion() || is_insertion();
  }

};


inline bool
operator<(Event2 const & a, Event2 const & b)
{
  // ordering is: insertion, deletion, snp
  long const order_a = a.is_deletion() + 2l * static_cast<long>(a.is_snp());
  long const order_b = b.is_deletion() + 2l * static_cast<long>(b.is_snp());

  return a.pos < b.pos ||
         (a.pos == b.pos && order_a < order_b) ||
         (a.pos == b.pos && order_a == order_b && a.ref < b.ref) ||
         (a.pos == b.pos && order_a == order_b && a.ref == b.ref && a.alt < b.alt);
}


inline bool
operator==(Event2 const & a, Event2 const & b)
{
  return a.pos == b.pos && a.ref == b.ref && a.alt == b.alt;
}


inline std::set<Event2>
get_edit_script(std::pair<std::string, std::string> const & s,
                bool const is_normalize,
                bool const is_trim_indel_on_ends)
{
  using Tedit_script = std::set<Event2>;
  assert(s.first.size() == s.second.size());

  Tedit_script edit_script;

  // reference string without gaps
  std::string ref;
  std::copy_if(s.first.begin(), s.first.end(), std::back_inserter(ref), [](char c){
      return c != '-';
    });

  // strings which are going to be used to create a new variant event
  std::vector<char> s1;
  std::vector<char> s2;

  // position of the new variant event
  long pos{0};
  long pos_q{0};

  auto are_s1_s2_empty =
    [&s1, &s2]() -> bool {
      return s1.empty() && s2.empty();
    };

  auto add_to_edit_script =
    [&pos, &pos_q, &s, &s1, &s2, &edit_script, &ref, is_normalize, is_trim_indel_on_ends]() -> void
    {
      long event_position = pos - static_cast<long>(s1.size());
      long event_position_q = pos_q - static_cast<long>(s2.size());
      assert(event_position >= 0);
      assert(event_position_q >= 0);

      if (is_trim_indel_on_ends && (s1.size() != s2.size() && event_position == 0))
      {
        s1.clear();
        s2.clear();
        return;
      }

      // Normalization is not possible if event_position is zero
      if (is_normalize && event_position > 0)
      {
        if (s1.empty())
        {
          // Normalize insertion
          long const s2_size = s2.size();
          assert(event_position - 1 < static_cast<long>(ref.size()));

          while (event_position > 0 && s2[s2_size - 1] == ref[event_position - 1])
          {
            s2.insert(s2.begin(), s2[s2_size - 1]); // Insert base in front
            s2.resize(s2_size); // Remove base in back to keep the same size
            --event_position; // Adjust event position accordingly
            --event_position_q;
          }
        }
        else if (s2.empty())
        {
          // Normalize deletion
          long const s1_size = s1.size();
          assert(event_position - 1 < static_cast<long>(ref.size()));

          while (event_position > 0 && s1[s1_size - 1] == ref[event_position - 1])
          {
            s1.insert(s1.begin(), s1[s1_size - 1]); // Insert base in front
            s1.resize(s1_size); // Remove base in back to keep the same size
            --event_position; // Adjust event position accordingly
            --event_position_q;
          }
        }
      }

      {
        // Create a new variant event
        Event2 new_edit(event_position, event_position_q,
                        std::string(s1.begin(), s1.end()),
                        std::string(s2.begin(), s2.end()));

        edit_script.insert(std::move(new_edit));
      }

      s1.clear();
      s2.clear();
    };

  for (long i = 0; i < static_cast<long>(s.first.size()); ++i)
  {
    // Both string cant have a gap
    assert(s.first[i] != '-' || s.second[i] != '-');

    if (s.first[i] == '-')
    {
      // Insertion
      if (are_s1_s2_empty() || (s2.size() > 0 && s1.size() == 0))
      {
        s2.push_back(s.second[i]); // Extend the insertion
      }
      else
      {
        add_to_edit_script(); // add previous variant
        s2.push_back(s.second[i]); // start new insertion
      }

      ++pos_q;
    }
    else if (s.second[i] == '-')
    {
      // Deletion
      if (are_s1_s2_empty() || (s1.size() > 0 && s2.size() == 0))
      {
        s1.push_back(s.first[i]); // Extend the deletion
      }
      else
      {
        add_to_edit_script(); // add previous variant
        s1.push_back(s.first[i]); // start new deletion
      }

      ++pos;
    }
    else if (s.first[i] != s.second[i] && s.first[i] != 'N' && s.second[i] != 'N')
    {
      // Mismatch
      if (!are_s1_s2_empty())
        add_to_edit_script();

      ++pos;
      ++pos_q;
      s1.push_back(s.first[i]);
      s2.push_back(s.second[i]);
      add_to_edit_script();
    }
    else
    {
      // Match
      if (!are_s1_s2_empty())
        add_to_edit_script();

      ++pos;
      ++pos_q;
    }
  }

  if (!are_s1_s2_empty() && !is_trim_indel_on_ends)
    add_to_edit_script();

  return edit_script;
}


} // namespace paw
