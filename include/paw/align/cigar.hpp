#pragma once

#include <cstdint>

namespace paw
{

enum struct CigarOperation
{
  MATCH = 0,
  INSERTION = 1,
  DELETION = 2,
  REF_SKIP = 3,
  SOFT_CLIP = 4,
  HARD_CLIP = 5,
  PAD = 6,
  EQUAL = 7,
  DIFFERENT = 8,
  BACK = 9,
  UNSET = 10
};

inline char cigar2char(CigarOperation op)
{
  return "MIDNSHP=XBU"[static_cast<int>(op)];
}

inline char inv_cigar2char(CigarOperation op)
{
  return "MDINSHP=XBU"[static_cast<int>(op)];
}

struct Cigar
{
  uint32_t count{0};
  CigarOperation operation{CigarOperation::UNSET};

  Cigar() = default;

  Cigar(uint32_t _count, CigarOperation op)
    : count(_count)
    , operation(op)
    {}
};

} // namespace paw
