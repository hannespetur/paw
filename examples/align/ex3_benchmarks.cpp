#include <paw/align.hpp>
#include <paw/parser.hpp>
#include <paw/internal/config.hpp>

#include <chrono> // std::chrono::high_resolution_clock
#include <set>
#include <string>
#include <vector>
#include <fstream> // std::ifstream

struct Test
{
  std::string seq1;
  std::string seq2;
  long expected_score;

  Test(std::string _seq1, std::string _seq2, long _expected_score)
    : seq1(_seq1)
    , seq2(_seq2)
    , expected_score(_expected_score)
  {}
};


int
main(int argc, char ** argv)
{
  std::vector<int> tests_to_run;

  try
  {
    paw::Parser parser(argc, argv);
    parser.set_name("Alignment example 3 - Running tests.");
    parser.parse_remaining_positional_arguments(tests_to_run,
                                                "list of tests to run...",
                                                "List of all tests to run (all if not specified)."
                                                );
    parser.finalize();
  }
  catch (const std::exception & e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  std::vector<Test> const tests =
  {
    {"GGG", "GGG", 6},//test 0
    {"GGGG", "GGG", 1},//test 1
    {"GGGGG", "GGG", 0},//test 2
    {"GGG", "GGGG", 1},//test 3
    {"GGG", "GGGGG", 0},//test 4
    {"AAA", "GGG", -6},//test 5
    {"CCCCCAAGGGGG", "CCCCCGGGGG", 14},//test 6
    {"TTTTTCCCCCAAGGGGGTTTTT", "TTTTTCCCCCGGGGGTTTTT", 34},//test 7
    {"AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 40},//test 8
    {"AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT", -40},//test 9
    {"AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
     "AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
     160},//test 10
    {"AAGTGTGTTAATTAATTAATGCTTGTAGGA", "GTTTATGTAGCTTATTCTATCCAAAGCAAT", -12},//test 11
    {"AAGTGTGTTAATTAATTAATGCTT", "TGTTAATTAATTAATGCTTGGCAAT", 19},//test 12
    {"GT", "GAT", -1},//test 13
    {"AAGACATCACGATG", "AAGACACCCCGCACG", 11}//test 14
  };

  // Run all tests if no specific tests are specified
  if (tests_to_run.size() == 0)
  {
    for (long i = 0; i < static_cast<long>(tests.size()); ++i)
      tests_to_run.insert(tests_to_run.end(), i);
  }

  std::string const current_archs = paw::get_current_arch();
  std::cerr << "Current archs are: " << current_archs << "\n";
  std::cerr << "== global_alignment ==\n";

  for (auto i : tests_to_run)
  {
    assert(i < static_cast<long>(tests.size()));
    auto const & test = tests[i];
    paw::AlignerOptions<int8_t> opts;
    auto score = paw::global_alignment(test.seq1, test.seq2, opts);

    if (score != test.expected_score)
    {
      std::cerr << "NOT OK. Score mismatch in test " << i
                << ". Got score " << score
                << " but I expected " << test.expected_score << std::endl;
    }
    else
    {
      std::cerr << "OK\n";
    }
  }

  std::cerr << "== global_alignment_score ==\n";

  for (auto i : tests_to_run)
  {
    auto const & test = tests[i];
    paw::AlignerOptions<int8_t> opts;
    auto score = paw::global_alignment_score(test.seq1, test.seq2, opts);

    if (score != test.expected_score)
    {
      std::cerr << "NOT OK. Score mismatch in test " << i
                << ". Got score " << score
                << " but I expected " << test.expected_score << std::endl;
    }
    else
    {
      std::cerr << "OK\n";
    }
  }
}
