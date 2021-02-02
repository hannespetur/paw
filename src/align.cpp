#define IMPLEMENT_PAW

#include <iomanip>
#include <iostream>
#include <vector>

#include <simdpp/simd.h>
#include <simdpp/dispatch/get_arch_linux_cpuinfo.h>

#include <paw/align/alignment_options.hpp>
#include <paw/align/alignment_results.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>
#include <paw/align/pairwise_alignment.hpp>
#include <paw/internal/config.hpp>


#define SIMDPP_USER_ARCH_INFO ::simdpp::get_arch_linux_cpuinfo()


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{


std::string
get_current_arch()
{
  simdpp::Arch current_arch = simdpp::this_compile_arch();
  std::ostringstream ss;
  const char * sep = "";

  if (current_arch == simdpp::Arch::NONE_NULL)
  {
    ss << sep << "NONE_NULL"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_SSE2))
  {
    ss << sep << "SSE2"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_SSE3))
  {
    ss << sep << "SSE3"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_SSE4_1))
  {
    ss << sep << "SSE4_1"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_AVX))
  {
    ss << sep << "AVX"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_AVX2))
  {
    ss << sep << "AVX2"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_FMA3))
  {
    ss << sep << "FMA3"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_FMA4))
  {
    ss << sep << "FMA4"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_XOP))
  {
    ss << sep << "XOP"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_AVX512F))
  {
    ss << sep << "AVX512F"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_AVX512BW))
  {
    ss << sep << "AVX512BW"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_AVX512DQ))
  {
    ss << sep << "AVX512DQ"; sep = ",";
  }

  if (static_cast<bool>(current_arch & simdpp::Arch::X86_AVX512VL))
  {
    ss << sep << "AVX512VL";
  }

  return ss.str();
}


} // namespace SIMDPP_ARCH_NAMESPACE


SIMDPP_MAKE_DISPATCHER_RET0(get_current_arch, std::string)

SIMDPP_MAKE_DISPATCHER((template <typename Tuint, typename Tseq>)
                         (< Tuint, Tseq >)
                         (void)
                         (set_query)
                         ((AlignmentOptions<Tuint>&)opt, (Tseq const &) seq)
                       )

SIMDPP_MAKE_DISPATCHER((template <typename Tseq, typename Tuint>)
                         (< Tseq, Tuint >)
                         (void)
                         (pairwise_alignment)
                         ((Tseq const &) x, (Tseq const &) y, (
                           AlignmentOptions<Tuint>&)z
                         )
                       )

SIMDPP_INSTANTIATE_DISPATCHER(
  (template void pairwise_alignment<std::string, uint8_t>(
     std::string const & s1, std::string const & s2,
     AlignmentOptions<uint8_t>&o)),
  (template void pairwise_alignment<std::string, uint16_t>(
     std::string const & s1, std::string const & s2,
     AlignmentOptions<uint16_t>&o))
  )

SIMDPP_INSTANTIATE_DISPATCHER(
  (template void pairwise_alignment<std::vector<char>, uint8_t>(
     std::vector<char> const & s1, std::vector<char> const & s2,
     AlignmentOptions<uint8_t>&o)),
  (template void pairwise_alignment<std::vector<char>, uint16_t>(
     std::vector<char> const & s1, std::vector<char> const & s2,
     AlignmentOptions<uint16_t>&o))
  )


} // namespace paw
