#define IMPLEMENT_PAW

#include <iomanip>
#include <iostream>
#include <vector>

#include <simdpp/simd.h>
#include <simdpp/dispatch/get_arch_linux_cpuinfo.h>

#include <paw/align/aligner_options.hpp>
#include <paw/align/global_alignment.hpp>
#include <paw/align/global_alignment_score.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>
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

SIMDPP_MAKE_DISPATCHER((template <typename Tseq, typename Tuint>)
                         (< Tseq, Tuint >)
                         (long)
                         (global_alignment)
                         ((Tseq const &) x, (Tseq const &)y, (
                           AlignerOptions<Tuint> const &)z
                         )
                       )

SIMDPP_MAKE_DISPATCHER((template <typename Tseq, typename Tuint>)
                         (< Tseq, Tuint >)
                         (long)
                         (global_alignment_score)
                         ((Tseq const &) x, (Tseq const &)y, (
                           AlignerOptions<Tuint> const &)z
                         )
                       )

SIMDPP_INSTANTIATE_DISPATCHER(
  (template long global_alignment<std::string, int8_t>(
     std::string const & s1, std::string const & s2,
     AlignerOptions<int8_t> const & o)),
  (template long global_alignment<std::string, int16_t>(
     std::string const & s1, std::string const & s2,
     AlignerOptions<int16_t> const & o)),
  (template long global_alignment_score<std::string, int8_t>(
     std::string const & s1, std::string const & s2,
     AlignerOptions<int8_t> const & o)),
  (template long global_alignment_score<std::string, int16_t>(
     std::string const & s1, std::string const & s2,
     AlignerOptions<int16_t> const & o))
  )


} // namespace paw
