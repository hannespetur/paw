#pragma once

#include <string> // std::string

#include <simdpp/simd.h>


namespace paw
{
//namespace SIMDPP_ARCH_NAMESPACE
//{

// An event is an single point mutation with two alleles, reference and alternative.
class Event2
{
public:
  long pos; // Position of the event.
  std::string ref; // Reference allele.
  std::string alt; // Alternative allele.

  bool is_snp() const;
  bool is_deletion() const;
  bool is_insertion() const;

};


bool operator<(Event2 const & a, Event2 const & b);
bool operator==(Event2 const & a, Event2 const & b);


//} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#if defined(IMPLEMENT_PAW)

#include <set>

namespace paw
{
//namespace SIMDPP_ARCH_NAMESPACE
//{

bool
Event2::is_snp() const
{
  return ref.size() == 1 && alt.size() == 1;
}


bool
Event2::is_deletion() const
{
  return alt.size() == 0;
}


bool
Event2::is_insertion() const
{
  return ref.size() == 0;
}


bool
operator<(Event2 const & a, Event2 const & b)
{
  return a.pos < b.pos ||
         (a.pos == b.pos && a.ref < b.ref) ||
         (a.pos == b.pos && a.ref == b.ref && a.alt < b.alt);
}


bool
operator==(Event2 const & a, Event2 const & b)
{
  return a.pos == b.pos && a.ref == b.ref && a.alt == b.alt;
}


std::set<Event2> inline
get_edit_script(std::pair<std::string, std::string> const & s, bool const is_normalize)
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
  long pos = 0;

  auto are_s1_s2_empty =
    [&s1, &s2]() -> bool {
      return s1.empty() && s2.empty();
    };

  auto add_to_edit_script =
    [&pos, &s, &s1, &s2, &edit_script, &ref, is_normalize]() -> void
    {
      long event_position = pos - static_cast<long>(s1.size());
      assert(event_position >= 0);

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
          }
        }
      }

      // Create a new variant event
      Event2 new_edit =
      {
        event_position,
        std::string(s1.begin(), s1.end()),
        std::string(s2.begin(), s2.end())
      };

      edit_script.insert(std::move(new_edit));

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
        s2.push_back(s.second[i]); // Extend the insertion
      else
        add_to_edit_script();
    }
    else if (s.second[i] == '-')
    {
      // Deletion
      if (are_s1_s2_empty() || (s1.size() > 0 && s2.size() == 0))
        s1.push_back(s.first[i]); // Extend the deletion
      else
        add_to_edit_script();

      ++pos;
    }
    else if (s.first[i] != s.second[i])
    {
      // Mismatch
      if (!are_s1_s2_empty())
        add_to_edit_script();

      ++pos;
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
    }
  }

  if (!are_s1_s2_empty())
    add_to_edit_script();

  return edit_script;
}


//} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#endif // IMPLEMENT_PAW
