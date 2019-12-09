#pragma once

#include <cstdint> // uint32_t
#include <map> // std::map<T1,T2>
#include <string> // std::string
#include <sstream> // std::stringstream
#include <vector> // std::vector<T>

#include <paw/align/event.hpp>


namespace paw
{
//namespace SIMDPP_ARCH_NAMESPACE
//{

class Variant
{
public:
  std::map<SIMDPP_ARCH_NAMESPACE::Event2, uint32_t /*allele index*/> event_to_allele;
  std::vector<uint16_t> calls;

public:
  using Tseqs = std::vector<std::string>;

  uint32_t pos;   // Position of the variant
  Tseqs seqs;   // Allele sequences of the variant. The first sequence is the reference allele.

  Variant(uint32_t _pos = -1, Tseqs const & _seqs = Tseqs());

  /// Read only methods
  bool has_sequences() const;
  bool is_deletion() const;
  bool is_insertion() const;
  bool is_snp() const;
  long get_max_del_reach() const;
  long get_max_reach() const;
  uint16_t get_call(long index) const;

  void print_events_to_alleles(std::ostream & ss) const;
  void print_seqs(std::ostream & ss) const;
  ///

  /// Class modifiers
  void add_base_to_front(std::string const & reference);
  void add_call(uint16_t const call);
  void clear();
  void add_event(SIMDPP_ARCH_NAMESPACE::Event2 const & e);
///

};


bool operator<(Variant const & a, Variant const & b);
bool operator==(Variant const & a, Variant const & b);

//} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#if defined(IMPLEMENT_PAW) || defined(__JETBRAINS_IDE__)


namespace paw
{
//namespace SIMDPP_ARCH_NAMESPACE
//{


Variant::Variant(uint32_t _pos, Variant::Tseqs const & _seqs)
  : pos(_pos)
  , seqs(_seqs)
{}


bool
Variant::has_sequences() const
{
  return seqs.size() > 1;
}


bool
Variant::is_deletion() const
{
  return has_sequences() && seqs[seqs.size() - 1].size() == 0;
}


bool
Variant::is_insertion() const
{
  return has_sequences() && seqs[0].size() == 0;
}


bool
Variant::is_snp() const
{
  return has_sequences() && seqs[0].size() == 1 && seqs[1].size() == 1;
}


uint16_t
Variant::get_call(long index) const
{
  assert(index < static_cast<long>(calls.size()));
  return calls[index];
}


long
Variant::get_max_del_reach() const
{
  if (seqs.size() < 2 || !is_deletion())
    return 0;

  return pos + seqs[0].size() - seqs[seqs.size() - 1].size();
}


long
Variant::get_max_reach() const
{
  if (seqs.size() < 2)
    return 0;

  if (is_deletion())
    return pos + seqs[0].size() - seqs[seqs.size() - 1].size();

  if (is_insertion())
    return pos;

  return pos + 1; //otherwise its a snp
}


void
Variant::print_events_to_alleles(std::ostream & ss) const
{
  for (auto const & event_allele : event_to_allele)
    ss << " " << event_allele.first.pos << " " << event_allele.second << "\n";
}


void
Variant::print_seqs(std::ostream & ss) const
{
  ss << "Seqs @ " << pos << ": ";

  for (long s = 0; s < static_cast<long>(seqs.size()); ++s)
  {
    if (s > 0)
      ss << ",";

    auto const & seq = seqs[s];

    if (seq.size() == 0)
      ss << "-";
    else
      ss << seq;
  }

  ss << "\n";
}


void
Variant::add_call(uint16_t const call)
{
  calls.push_back(call);
}


void
Variant::add_base_to_front(std::string const & reference)
{
  --pos;

  for (auto & seq : seqs)
  {
    if (seq != "*")
      seq = reference[pos] + seq;
  }
}


void
Variant::clear()
{
  pos = -1;
  seqs.clear();
}


void
Variant::add_event(SIMDPP_ARCH_NAMESPACE::Event2 const & e)
{
  event_to_allele[e] = this->seqs.size();
  this->seqs.push_back(e.alt);
}


bool
operator<(Variant const & a, Variant const & b)
{
  return a.pos < b.pos || (a.pos == b.pos && a.seqs < b.seqs);
}


bool
operator==(Variant const & a, Variant const & b)
{
  return a.pos == b.pos && a.seqs == b.seqs;
}


//} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw

#endif // IMPLEMENT_PAW
