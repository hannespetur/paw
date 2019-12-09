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
  std::size_t pos;   // Position of the event.
  std::string ref;   // Reference allele.
  std::string alt;   // Alternative allele.

  bool is_snp() const;
  bool is_deletion() const;
  bool is_insertion() const;

};


bool operator<(Event2 const & a, Event2 const & b);
bool operator==(Event2 const & a, Event2 const & b);


//} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#if defined(IMPLEMENT_PAW)


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
get_edit_script(std::pair<std::string, std::string> const & s)
{
  using Tedit_script = std::set<Event2>;

  Tedit_script edit_script;
  std::vector<char> s1;
  std::vector<char> s2;
  std::size_t pos = 0;

  auto are_sequences_empty = [&](){
                               return s1.size() == 0 && s2.size() == 0;
                             };

  auto add_to_edit_script = [&]()
                            {
                              Event2 new_edit =
                              {
                                pos - s1.size(),
                                std::string(s1.begin(), s1.end()),
                                std::string(s2.begin(), s2.end())
                              };

                              edit_script.insert(std::move(new_edit));

                              s1.clear();
                              s2.clear();
                            };


  for (std::size_t i = 0; i < s.first.size(); ++i)
  {
    // Both string cant have a gap
    assert(s.first[i] != '-' || s.second[i] != '-');

    if (s.first[i] == '-')
    {
      // Insertion
      if (are_sequences_empty() || (s2.size() > 0 && s1.size() == 0))
        s2.push_back(s.second[i]); // Extend the insertion
      else
        add_to_edit_script();
    }
    else if (s.second[i] == '-')
    {
      // Deletion
      if (are_sequences_empty() || (s1.size() > 0 && s2.size() == 0))
        s1.push_back(s.first[i]); // Extend the deletion
      else
        add_to_edit_script();

      ++pos;
    }
    else if (s.first[i] != s.second[i])
    {
      // Mismatch
      if (!are_sequences_empty())
        add_to_edit_script();

      ++pos;
      s1.push_back(s.first[i]);
      s2.push_back(s.second[i]);
      add_to_edit_script();
    }
    else
    {
      // Match
      if (!are_sequences_empty())
        add_to_edit_script();

      ++pos;
    }
  }

  if (!are_sequences_empty())
    add_to_edit_script();

  return edit_script;
}


//} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#endif // IMPLEMENT_PAW
