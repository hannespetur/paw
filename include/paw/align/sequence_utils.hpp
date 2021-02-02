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


} // namespace paw
