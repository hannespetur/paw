#pragma once

#include <sstream>
#include <string>


namespace paw
{

inline std::string
reverse_complement(std::string const & sequence)
{
  std::ostringstream r;

  for (auto it = sequence.rbegin(); it != sequence.rend(); ++it)
  {
    switch (*it)
    {
    case 'A':
    case 'a':
      r << 'T'; break;

    case 'C':
    case 'c':
      r << 'G'; break;

    case 'G':
    case 'g':
      r << 'C'; break;

    case 'T':
    case 't':
      r << 'A'; break;

    case 'N':
    case 'n':
      r << 'N'; break;

    default:
      std::cerr << "Unexpected base: " << *it << std::endl;
      std::exit(1);

    }
  }

  return r.str();
}


inline void
remove_common_prefix(long & pos,
                     std::vector<std::string> & seqs,
                     bool keep_one_match = false)
{
  if (seqs.size() <= 1 || seqs[0].size() <= 1)
    return;

  auto & ref = seqs[0];
  long p{0}; // prefix length

  while ((static_cast<long>(ref.size()) - p) > 1)
  {
    char ref_base = ref[p];

    if (std::any_of(std::next(seqs.begin()),
                    seqs.end(),
                    [&](std::string const & alt)
      {
        return (static_cast<long>(alt.size()) - p) <= 1 ||
        alt[p] != ref_base ||
        (keep_one_match && alt[p + 1] != ref[p + 1]);
      }))
    {
      break;
    }

    ++p;
  }

  if (p > 0)
  {
    pos += p;

    for (auto & seq : seqs)
      seq.erase(seq.begin(), seq.begin() + p);
  }
}


inline void
remove_common_suffix(std::vector<std::string> & seqs)
{
  if (seqs.size() <= 1 || seqs[0].size() <= 1)
    return;

  auto & ref = seqs[0];
  char const last_ref = ref.back();

  while (ref.size() > 1)
  {
    for (long a{1}; a < static_cast<long>(seqs.size()); ++a)
    {
      auto const & alt = seqs[a];
      assert(alt.size() > 0);

      if (alt.size() <= 1 || alt[alt.size() - 1] != last_ref)
        return;
    }

    // Remove one base from back of reference
    for (auto & seq : seqs)
      seq.pop_back();
  }
}


} // namespace paw
