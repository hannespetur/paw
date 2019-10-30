#pragma once

#include <string> // std::string

#include <simdpp/simd.h>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

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


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#if defined(IMPLEMENT_PAW)


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

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


} // namespace paw
} // namespace SIMDPP_ARCH_NAMESPACE

#endif // IMPLEMENT_PAW
