#include <iostream>
#include <vector>

#include <simdpp/simd.h>
#include <simdpp/dispatch/get_arch_linux_cpuinfo.h>
#include <simdpp/dispatch/get_arch_string_list.h>

#include <paw/align/libsimdpp_align.hpp>
//#include <paw/align/libsimdpp_row.hpp>


#define SIMDPP_USER_ARCH_INFO ::simdpp::get_arch_linux_cpuinfo()


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

namespace simd = simdpp;


void
align(std::string const & seq1, std::string const & seq2)
{
  std::cerr << "Aligning " << seq1 << " to " << seq2 << "\n";

  Align<std::string::const_iterator> align(seq1.cbegin(), seq2.cend());
  align.calculate_DNA_W_profile();
}


std::string
get_current_arch()
{
  simd::Arch current_arch = simd::this_compile_arch();
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
    ss << sep << "AVX512VL"; sep = ",";
  }

  return ss.str();

  /*
  simd::int16<8> test3 = simd::make_zero();
  simd::int16<8> test4 = simd::make_int(1, 0, 0, 0, 0, 0, 0, 42);
  auto test5 = add(test3, test4);

  //simd::int16v vec = test3.vec(0);
  constexpr int test = test3.vec_length / test3.length;
  std::cout << test << " "
            << test3.vec_length << " "
            << test3.length << " "
            << test3.base_length << " "
            << test3.vec(0).vec_length << " "
            << test3.vec(0).length << " "
            << test3.vec(0).base_length << " "
            << std::endl;

  std::vector<int16_t> result(8);
  simdpp::store(&result[0], test5);

  for (auto const & r : result)
    std::cout << r << " ";

  std::cout << std::endl;
  */
}


} // namespace SIMDPP_ARCH_NAMESPACE

SIMDPP_MAKE_DISPATCHER_RET0(get_current_arch, std::string)
SIMDPP_MAKE_DISPATCHER_RET2(align, void, std::string const &, std::string const &)

} // namespace paw
