#pragma once

#include <cstdint>


namespace paw
{

enum CigarOperation
{
  MATCH = 0,
  INSERTION,
  DELETION
};


struct Cigar
{
  std::size_t count;
  CigarOperation operation;
};

} // namespace paw
