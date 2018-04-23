#define IMPLEMENT_PAW

#include <iostream>
#include <vector>

#include <simdpp/simd.h>
#include <simdpp/dispatch/get_arch_linux_cpuinfo.h>

#include <paw/align/libsimdpp_align.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/internal/config.hpp>


#define SIMDPP_USER_ARCH_INFO ::simdpp::get_arch_linux_cpuinfo()


namespace paw
{

namespace SIMDPP_ARCH_NAMESPACE
{

long
align(std::string const & seq1, std::string const & seq2)
{
  std::string const nspace(STR(SIMDPP_ARCH_NAMESPACE));
  std::cerr << "Namespace = " << nspace << "\n";
  Align<std::string::const_iterator> align(seq1.cbegin(), seq1.cend());
  align.calculate_DNA_W_profile();
  auto score = align.align(seq2.cbegin(), seq2.cend());
  auto aligned_strings = align.get_aligned_strings();

  for (long i = 0; i < std::min(1000l, static_cast<long>(aligned_strings.first.size())); i += 140)
  {
    std::cout << aligned_strings.first.substr(i, 140) << "\n"
              << aligned_strings.second.substr(i, 140) << "\n\n";
  }

  return score;
}


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

  /*
  simdpp::int16<8> test3 = simdpp::make_zero();
  simdpp::int16<8> test4 = simdpp::make_int(1, 0, 0, 0, 0, 0, 0, 42);
  auto test5 = add(test3, test4);

  //simdpp::int16v vec = test3.vec(0);
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
SIMDPP_MAKE_DISPATCHER_RET2(align, long, std::string const &, std::string const &)

} // namespace paw
