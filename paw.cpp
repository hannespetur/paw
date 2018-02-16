#define IMPLEMENT_PAW
#include "paw.hpp"

template class paw::Aligner<uint8_t, std::string::const_iterator>;
template class paw::Aligner<uint16_t, std::string::const_iterator>;
template class paw::Aligner<uint32_t, std::string::const_iterator>;
//template class paw::Aligner<std::string::iterator, uint64_t>;

template struct paw::AlignerOptions<uint8_t>; //::get_score_T(char const b) const;
template struct paw::AlignerOptions<uint16_t>; //::get_score_T(char const b) const;
template struct paw::AlignerOptions<uint32_t>; //::get_score_T(char const b) const;
//template struct paw::AlignerOptions<uint64_t>; //::get_score_T(char const b) const;

template struct paw::Backtracker<uint8_t>;
template struct paw::Backtracker<uint16_t>;
template struct paw::Backtracker<uint32_t>;
//template struct paw::Backtracker<uint64_t>;

template struct paw::Row<uint8_t>;
template struct paw::Row<uint16_t>;
template struct paw::Row<uint32_t>;
