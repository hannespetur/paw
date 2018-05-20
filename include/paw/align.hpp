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

template <typename Tseq, typename Tuint>
long
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignerOptions<Tuint> const & _opt =
                   AlignerOptions<Tuint>(true /*default options*/)
                 );

template <typename Tseq, typename Tuint>
long
global_alignment_score(Tseq const & seq1,
                       Tseq const & seq2,
                       AlignerOptions<Tuint> const & _opt =
                         AlignerOptions<Tuint>(true /*default options*/)
                       );


// See: https://github.com/p12tic/libsimdpp/blob/3ab0be0e5aa0773f152d7d759400173d64253534/simdpp/detail/insn_id.h

namespace arch_null
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
long
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignerOptions<Tuint> const & _opt =
                   AlignerOptions<Tuint>(true /*default options*/)
                 );

template <typename Tseq, typename Tuint>
long
global_alignment_score(Tseq const & seq1,
                       Tseq const & seq2,
                       AlignerOptions<Tuint> const & _opt =
                         AlignerOptions<Tuint>(true /*default options*/)
                       );


}

namespace arch_sse2
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
long
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignerOptions<Tuint> const & _opt =
                   AlignerOptions<Tuint>(true /*default options*/)
                 );

template <typename Tseq, typename Tuint>
long
global_alignment_score(Tseq const & seq1,
                       Tseq const & seq2,
                       AlignerOptions<Tuint> const & _opt =
                         AlignerOptions<Tuint>(true /*default options*/)
                       );


}

namespace arch_sse3
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
long
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignerOptions<Tuint> const & _opt = AlignerOptions<Tuint>(true /*default options*/)
                 );

template <typename Tseq, typename Tuint>
long
global_alignment_score(Tseq const & seq1,
                       Tseq const & seq2,
                       AlignerOptions<Tuint> const & _opt =
                         AlignerOptions<Tuint>(true /*default options*/)
                       );


}

namespace arch_sse4p1
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
long
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignerOptions<Tuint> const & _opt =
                   AlignerOptions<Tuint>(true /*default options*/)
                 );

template <typename Tseq, typename Tuint>
long
global_alignment_score(Tseq const & seq1,
                       Tseq const & seq2,
                       AlignerOptions<Tuint> const & _opt =
                         AlignerOptions<Tuint>(true /*default options*/)
                       );


}

namespace arch_popcnt_avx
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
long
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignerOptions<Tuint> const & _opt =
                   AlignerOptions<Tuint>(true /*default options*/)
                 );

template <typename Tseq, typename Tuint>
long
global_alignment_score(Tseq const & seq1,
                       Tseq const & seq2,
                       AlignerOptions<Tuint> const & _opt =
                         AlignerOptions<Tuint>(true /*default options*/)
                       );


}


namespace arch_popcnt_avx2
{

std::string get_current_arch();

template <typename Tseq, typename Tuint>
long
global_alignment(Tseq const & seq1,
                 Tseq const & seq2,
                 AlignerOptions<Tuint> const & _opt =
                   AlignerOptions<Tuint>(true /*default options*/)
                 );

template <typename Tseq, typename Tuint>
long
global_alignment_score(Tseq const & seq1,
                       Tseq const & seq2,
                       AlignerOptions<Tuint> const & _opt =
                         AlignerOptions<Tuint>(true /*default options*/)
                       );

}

} // namespace paw
