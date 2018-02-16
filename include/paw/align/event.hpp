#pragma once

#include <string> // std::string


namespace paw
{

// An event is an single point mutation with two alleles, reference and alternative.
class Event
{
public:
  std::size_t pos;   // Position of the event.
  std::string ref;   // Reference allele.
  std::string alt;   // Alternative allele.

  bool is_snp() const;
  bool is_deletion() const;
  bool is_insertion() const;

};


bool operator<(Event const & a, Event const & b);
bool operator==(Event const & a, Event const & b);


} // namespace paw


#ifdef IMPLEMENT_PAW


namespace paw
{


bool
Event::is_snp() const
{
  return ref.size() == 1 && alt.size() == 1;
}


bool
Event::is_deletion() const
{
  return alt.size() == 0;
}


bool
Event::is_insertion() const
{
  return ref.size() == 0;
}


bool
operator<(Event const & a, Event const & b)
{
  return a.pos < b.pos ||
         (a.pos == b.pos && a.ref < b.ref) ||
         (a.pos == b.pos && a.ref == b.ref && a.alt < b.alt);
}


bool
operator==(Event const & a, Event const & b)
{
  return a.pos == b.pos && a.ref == b.ref && a.alt == b.alt;
}


} // namespace paw


#endif // IMPLEMENT_PAW
