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

inline CigarOperation inv_cigar(CigarOperation op)
{
  switch(op)
  {
  case CigarOperation::DELETION:
    return CigarOperation::INSERTION;
  case CigarOperation::INSERTION:
    return CigarOperation::DELETION;
  default:
    return op;
  }
}

inline bool advances_ref(CigarOperation op)
{
  switch(op)
  {
  case CigarOperation::MATCH:
  case CigarOperation::DELETION:
  case CigarOperation::REF_SKIP:
  case CigarOperation::EQUAL:
  case CigarOperation::DIFFERENT:
    return true;
  default:
    return false;
  }
}

inline bool advances_query(CigarOperation op)
{
  switch(op)
  {
  case CigarOperation::MATCH:
  case CigarOperation::INSERTION:
  case CigarOperation::SOFT_CLIP:
  case CigarOperation::EQUAL:
  case CigarOperation::DIFFERENT:
    return true;
  default:
    return false;
  }
}

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
