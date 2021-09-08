#pragma once

#include <cstdint>


namespace paw
{

enum CigarOperation
{
  MATCH = 0,
  INSERTION,
  DELETION,
  UNSET
};


struct Cigar
{
  uint32_t count{0};
  CigarOperation operation{UNSET};
};

} // namespace paw
