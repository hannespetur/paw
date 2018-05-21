#include <paw/align.hpp>
#include <paw/internal/config.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <chrono> // std::chrono::high_resolution_clock
#include <string>
#include <vector>
#include <fstream> // std::ifstream


namespace io = boost::iostreams;

namespace
{

inline std::string
get_sequence_from_fa(std::string const & fn, bool gzip = false)
{
  std::ifstream file(fn, std::ios_base::in | std::ios_base::binary);

  io::filtering_istream in;

  if (gzip)
    in.push(io::gzip_decompressor());

  in.push(file);

  std::stringstream ss;
  std::string str;
  std::getline(in, str); // Throw away header

  for (; std::getline(in, str); )
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
  std::string database = get_sequence_from_fa(src_dir + "/test/data/MT-human.fa", false);
  //database = get_sequence_from_fa(PROJECT_SOURCE_DIR + "/test/data/t2.fa.gz", true);

  std::string query = get_sequence_from_fa(src_dir + "/test/data/MT-orang.fa", false);

  // Change to uppercase
  std::transform(database.begin(), database.end(), database.begin(), ::toupper);
  std::transform(query.begin(), query.end(), query.begin(), ::toupper);

  using Ttime = std::chrono::high_resolution_clock;
  using Tduration = std::chrono::duration<double, std::milli>;
  paw::AlignerOptions<int16_t> opts(true /*default options*/);

  //database = database.substr(200, 25);
  //query = query.substr(0, 25);

  auto t0 = Ttime::now();
  auto score = paw::arch_null::global_alignment(database, query, opts);
  auto t1 = Ttime::now();
  std::cout << "null " << Tduration(t1 - t0).count() << " ms\n";
  std::cout << "score = " << score << "\n";

  t0 = Ttime::now();
  score = paw::arch_sse2::global_alignment(database, query, opts);
  t1 = Ttime::now();
  std::cout << "sse2 " << Tduration(t1 - t0).count() << " ms\n";
  std::cout << "score = " << score << "\n";

  t0 = Ttime::now();
  score = paw::arch_sse3::global_alignment(database, query, opts);
  t1 = Ttime::now();
  std::cout << "sse3 " << Tduration(t1 - t0).count() << " ms\n";
  std::cout << "score = " << score << "\n";

  t0 = Ttime::now();
  score = paw::arch_sse4p1::global_alignment(database, query, opts);
  t1 = Ttime::now();
  std::cout << "sse4_1 " << Tduration(t1 - t0).count() << " ms\n";
  std::cout << "score = " << score << "\n";

  t0 = Ttime::now();
  score = paw::global_alignment(database, query, opts);
  t1 = Ttime::now();
  std::cout << "best " << Tduration(t1 - t0).count() << " ms\n";
  std::cout << "score = " << score << "\n";
}
