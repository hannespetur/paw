#pragma once

#include <string>

#include <paw/align/libsimdpp_align.hpp>

//#include <paw/align/boost_simd_align.hpp>
//#include <paw/align/event.hpp>
//#include <paw/align/fasta.hpp>
//#include <paw/align/row.hpp>
//#include <paw/align/seq_it.hpp>
//#include <paw/align/skyr.hpp>
//#include <paw/align/variant.hpp>
//#include <paw/align/vcf.hpp>

namespace paw
{

std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);

// See: https://github.com/p12tic/libsimdpp/blob/3ab0be0e5aa0773f152d7d759400173d64253534/simdpp/detail/insn_id.h

namespace arch_null
{
std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);
}

namespace arch_sse2
{
std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);
}

namespace arch_sse3
{
std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);
}

namespace arch_sse4p1
{
std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);
}

namespace arch_avx_popcnt
{
std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);
}

namespace arch_avx_fm4
{
std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);
}

namespace arch_avx2_fm3
{
std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);
}
} // namespace paw
