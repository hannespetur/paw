#pragma once

#include <string>

#include <paw/align/alignment_results.hpp>
#include <paw/align/alignment_options.hpp>
#include <paw/align/libsimdpp_align.hpp>


namespace paw
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );


/// See namespaces at:
/// https://github.com/p12tic/libsimdpp/blob/3ab0be0e5aa0773f152d7d759400173d64253534/simdpp/detail/insn_id.h

namespace arch_null
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );

}

namespace arch_sse2
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );

}

namespace arch_sse3
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );

}

namespace arch_sse4p1
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );

}

namespace arch_sse4p1_popcnt
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );

}

namespace arch_popcnt_avx
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );

}


namespace arch_popcnt_avx2
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
void
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignmentOptions<Tuint> & opts
                 );

}

} // namespace paw
