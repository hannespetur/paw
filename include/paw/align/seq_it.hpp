#pragma once

#include <string> // std::string
#include <iostream> // std::ostream
#include <vector> // std::vector<T>


namespace paw
{

struct SeqIt // begin and end iterators for each sequence
{
  using Tit = std::string::const_iterator;
  Tit begin;
  Tit end;
  Tit it;

  SeqIt(Tit _begin, Tit _end)
    : begin(_begin)
    , end(_end)
    , it(_begin)
  {}


  SeqIt &
  operator++()
  {
    if (it != end)
      ++it;

    return *this;
  }


  const char *
  get() const
  {
    if (!is_at_end())
      return &(*it);

    return nullptr;
  }


  friend std::ostream &
  operator<<(std::ostream & output, SeqIt const & a)
  {
    if (!a.is_at_end())
      output << *a.it;

    return output;
  }


  bool inline
  is_at_end() const
  {
    return it == end;
  }


};


std::vector<SeqIt> inline
make_seq_its(std::vector<std::string> const & seqs)
{
  std::vector<SeqIt> seq_its;
  seq_its.reserve(seqs.size());

  for (auto const & seq : seqs)
    seq_its.push_back({seq.cbegin(), seq.cend()});

  return seq_its;
}


} // namespace paw
