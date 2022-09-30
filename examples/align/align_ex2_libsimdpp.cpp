#include <paw/align.hpp>
#include <paw/internal/config.hpp>

#include <chrono> // std::chrono::high_resolution_clock
#include <string>
#include <vector>
#include <fstream> // std::ifstream


namespace
{

std::string
get_sequence_from_fa(std::string const & fn)
{
  std::ifstream in(fn, std::ios_base::in | std::ios_base::binary);

  std::stringstream ss;
  std::string str;
  std::getline(in, str); // Throw away header

  for (; std::getline(in, str);)
  {
    if (str[0] == '>')
      break;

    ss << str;
  }

  return ss.str();
}


} // anon namespace


int
main(int, char **)
{
  std::string const current_archs = paw::get_current_arch();
  std::cerr << "Current archs are: " << current_archs << "\n";

  std::string const src_dir(STR(PAW_SOURCE_DIR));
  std::string database = get_sequence_from_fa(src_dir + "/test/data/MT-human.fa");
  //database = get_sequence_from_fa(PROJECT_SOURCE_DIR + "/test/data/t2.fa.gz", true);

  std::string query = get_sequence_from_fa(src_dir + "/test/data/MT-orang.fa");

  // Change to uppercase
  std::transform(database.begin(), database.end(), database.begin(), ::toupper);
  std::transform(query.begin(), query.end(), query.begin(), ::toupper);

  using Ttime = std::chrono::high_resolution_clock;
  using Tduration = std::chrono::duration<double, std::milli>;
  paw::AlignmentOptions<uint16_t> opts;

  //database = database.substr(200, 25);
  //query = query.substr(0, 25);

  auto t0 = Ttime::now();
  //paw::arch_null::pairwise_alignment(database, query, opts);
  auto t1 = Ttime::now();
  //std::cout << "null " << Tduration(t1 - t0).count() << " ms\n";
  //std::cout << "score = " << opts.get_alignment_results()->score << "\n";

  //t0 = Ttime::now();
  //opts = paw::AlignmentOptions<uint16_t>();
  //paw::arch_popcnt_avx2::pairwise_alignment(database, query, opts);
  //t1 = Ttime::now();
  //std::cout << "avx2 " << Tduration(t1 - t0).count() << " ms\n";
  //std::cout << "score = " << opts.get_alignment_results()->score << "\n";


  t0 = Ttime::now();
  opts = paw::AlignmentOptions<uint16_t>();

#ifdef __x86_64__
    paw::arch_sse2::pairwise_alignment(database, query, opts);
    std::string arch = "sse2 ";
#elif __aarch64__
    paw::arch_neon::pairwise_alignment(database, query, opts);
    std::string arch = "neon ";
#endif

  t1 = Ttime::now();
  std::cout << arch << Tduration(t1 - t0).count() << " ms\n";
  std::cout << "score = " << opts.get_alignment_results()->score << "\n";


  t0 = Ttime::now();
  opts = paw::AlignmentOptions<uint16_t>();
  paw::pairwise_alignment(database, query, opts);
  t1 = Ttime::now();
  std::cout << "best " << Tduration(t1 - t0).count() << " ms\n";
  std::cout << "score = " << opts.get_alignment_results()->score << "\n";
}
