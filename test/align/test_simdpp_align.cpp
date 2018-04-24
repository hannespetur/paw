#include "../include/catch.hpp"

#include <simdpp/simd.h>

#include <paw/align/libsimdpp_align.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <chrono> // std::chrono::high_resolution_clock
#include <string>
#include <vector>
#include <fstream> // std::ifstream


/*
TEST_CASE("W profile")
{
  using namespace paw::SIMDPP_ARCH_NAMESPACE;

  std::string database = "GCAG";
  std::string query =    "GAAG";

  using Tuint = uint16_t;
  paw::AlignerOptions<Tuint> opt;
  opt.match = 2;
  opt.mismatch = 2;
  opt.gap_open = 5;
  opt.gap_extend = 1;
  opt.backtracking = true;
  opt.top_row_free = false;
  opt.bottom_row_free = false;
  opt.top_row_gap_open_free = false;
  opt.bottom_row_gap_open_free = false;
  opt.left_column_gap_open_free = false;
  opt.right_column_gap_open_free = false;

  Align<std::string::const_iterator> align(database.cbegin(), database.cend(), opt);
  align.calculate_DNA_W_profile();

  SECTION("W profile contains 4 values for each DNAbase")
  {
    auto const & W_profile = align.get_W_profile();

    REQUIRE(W_profile.size() == 4);
    REQUIRE(W_profile[0].vectors.size() == 1);
    REQUIRE(W_profile[1].vectors.size() == 1);
    REQUIRE(W_profile[2].vectors.size() == 1);
    REQUIRE(W_profile[3].vectors.size() == 1);

    std::vector<Tuint, simdpp::aligned_allocator<Tuint, sizeof(Tuint)> > elements;
    elements.resize(16 / sizeof(Tuint) * W_profile[0].vectors[0].vec_length, 100);

    simdpp::store(&elements[0], W_profile[0].vectors[0]);
    REQUIRE(elements[0] == 0);
    REQUIRE(elements[1] == 0);
    REQUIRE(elements[2] == 4);
    REQUIRE(elements[3] == 0);

    simdpp::store(&elements[0], W_profile[1].vectors[0]);
    REQUIRE(elements[0] == 0);
    REQUIRE(elements[1] == 4);
    REQUIRE(elements[2] == 0);
    REQUIRE(elements[3] == 0);

    simdpp::store(&elements[0], W_profile[2].vectors[0]);
    REQUIRE(elements[0] == 4);
    REQUIRE(elements[1] == 0);
    REQUIRE(elements[2] == 0);
    REQUIRE(elements[3] == 4);

    simdpp::store(&elements[0], W_profile[3].vectors[0]);
    REQUIRE(elements[0] == 0);
    REQUIRE(elements[1] == 0);
    REQUIRE(elements[2] == 0);
    REQUIRE(elements[3] == 0);
  }

  align.align(query.cbegin(), query.cend());
}
*/


/*
TEST_CASE("simdpp alignment black box tests")
{
  using namespace paw;

  using Ttime = std::chrono::high_resolution_clock;
  using Tduration = std::chrono::duration<double, std::milli>;

  std::string database = "GGGG";
  std::string query =    "GGG";
  //std::swap(database, query);

  // Black box test with FASTA sequences runs too slowly in debug mode
#ifdef NDEBUG
  std::string const src_dir(STR(PAW_SOURCE_DIR));
  database = get_sequence_from_fa(src_dir + "/test/data/MT-human.fa", false);
  //database = get_sequence_from_fa(PROJECT_SOURCE_DIR + "/test/data/t2.fa.gz", true);

  query = get_sequence_from_fa(src_dir + "/test/data/MT-orang.fa", false);
  //query = get_sequence_from_fa(PROJECT_SOURCE_DIR + "/test/data/q2.fa.gz", true);

  // Change to uppercase
  std::transform(database.begin(), database.end(), database.begin(), ::toupper);
  std::transform(query.begin(), query.end(), query.begin(), ::toupper);
#endif // NDEBUG

  //boost_simd_align<int8_t>(query, database);
  auto t0 = Ttime::now();

  {
    paw::AlignerOptions<uint8_t> opt;
    opt.match = 2;
    opt.mismatch = 2;
    opt.gap_open = 5;
    opt.gap_extend = 1;
    opt.backtracking = false;
    opt.top_row_free = false;
    opt.bottom_row_free = false;
    opt.top_row_gap_open_free = false;
    opt.bottom_row_gap_open_free = false;
    opt.left_column_gap_open_free = false;
    opt.right_column_gap_open_free = false;
    paw::Align<std::string::const_iterator> aligner(database.cbegin(), database.cend(), opt);
    auto score = aligner.align(query.cbegin(), query.cend());
    std::cout << "score = " << score << "\n";
  }

  auto t1 = Ttime::now();
  //std::cout << "int8_t  " << Tduration(t1 - t0).count() << " ms\n";

  for (std::size_t i = 0; i < 1; ++i)
  {
    using Tint = uint16_t;
    t0 = Ttime::now();
    AlignerOptions<Tint> opt;
    opt.match = 2;
    opt.mismatch = 2;
    opt.gap_open = 5;
    opt.gap_extend = 1;
    opt.backtracking = true;
    opt.top_row_free = false;
    opt.bottom_row_free = false;
    opt.top_row_gap_open_free = false;
    opt.bottom_row_gap_open_free = false;
    opt.left_column_gap_open_free = false;
    opt.right_column_gap_open_free = false;
    paw::Align<std::string::const_iterator> align(database.cbegin(), database.cend(), opt);
    align.calculate_DNA_W_profile();
    auto score = align.align(query.cbegin(), query.cend());

    t1 = Ttime::now();
    std::cout << "score = " << score << "\n";
    //std::swap(database, query);
    std::cout << "int16_t " << Tduration(t1 - t0).count() << " ms\n";

    if (opt.backtracking)
    {
      auto aligned_strings = align.get_aligned_strings();

      for (std::size_t i = 0; i < std::min(1000ul, aligned_strings.first.size()); i += 90)
      {
        std::cout << aligned_strings.first.substr(i, 90) << "\n"
                  << aligned_strings.second.substr(i, 90) << "\n\n";
      }

      //auto edit_script = get_edit_script(aligned_strings);
      //
      //for (auto const & e : edit_script)
      //{
      //  std::cout << e.pos << " "
      //            << (e.ref.size() > 0 ? std::string(e.ref.begin(), e.ref.end()) : "-") << " "
      //            << (e.alt.size() > 0 ? std::string(e.alt.begin(), e.alt.end()) : "-") << "\n";
      //}

      auto cigar = align.get_cigar(aligned_strings);
      std::size_t M = 0;
      std::size_t I = 0;
      std::size_t D = 0;

      for (auto const & c : cigar)
      {
        std::cout << c.count;

        switch (c.operation)
        {
        case MATCH:
          std::cout << "M ";
          M += c.count;
          break;

        case INSERTION:
          std::cout << "I ";
          I += c.count;
          break;

        case DELETION:
          std::cout << "D ";
          D += c.count;
          break;

        default:
          std::cout << "WHAT\n";
          break;
        }
      }

      std::cout << "\n";
      std::cout << "M,I,D = " << M << "," << I << "," << D << "\n";
    }
  }
}
*/
