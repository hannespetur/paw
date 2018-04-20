#pragma once

#include <paw/align/backtracker.hpp>
#include <paw/align/boost_simd_align.hpp>
#include <paw/align/event.hpp>
#include <paw/align/fasta.hpp>
#include <paw/align/row.hpp>
#include <paw/align/seq_it.hpp>
#include <paw/align/skyr.hpp>
#include <paw/align/variant.hpp>
#include <paw/align/vcf.hpp>

namespace paw
{

std::string get_current_arch();
long align(std::string const & seq1, std::string const & seq2);

} // namespace paw
