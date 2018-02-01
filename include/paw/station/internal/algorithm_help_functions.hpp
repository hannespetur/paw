#pragma once


namespace paw
{
namespace stations_internal
{

/** Returns a number with a single bit set, the bit set is the most significant bit of the input. If the input is 0, then 1 is returned. Examples:
 * highest_ordered_bit(0) = 1
 * highest_ordered_bit(1) = 1
 * highest_ordered_bit(2) = 2
 * highest_ordered_bit(3) = 2
 * highest_ordered_bit(4) = 4
 * highest_ordered_bit(7) = 4
 * highest_ordered_bit(9) = 8
 */
template <typename T>
T inline
highest_ordered_bit(T num)
{
  T ret = 1;

  while (num >>= 1)
    ret <<= 1;

  return ret;
}


} // namespace stations_internal
} // namespace paw
