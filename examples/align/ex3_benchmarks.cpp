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
  std::string seq1{};
  std::string seq2{};
  long expected_score = 0;
  uint8_t match = 2;
  uint8_t mismatch = 2;
  uint8_t gap_open = 5;
  uint8_t gap_extend = 1;

  Test() = default;

  Test(std::string _seq1,
       std::string _seq2,
       long _expected_score = 0,
       uint8_t score_match = 2,
       uint8_t score_mismatch = 2,
       uint8_t score_gap_open = 5,
       uint8_t score_gap_extend = 1
       )
    : seq1(_seq1)
    , seq2(_seq2)
    , expected_score(_expected_score)
    , match(score_match)
    , mismatch(score_mismatch)
    , gap_open(score_gap_open)
    , gap_extend(score_gap_extend)
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
    {"GGG", "GGG", 6, 2 /*match*/, 2 /*mismatch*/, 10 /*gap_open*/, 1 /*gap extend*/}, //test 0
    {"GGGG", "GGG", 1}, //test 1
    {"GGGGG", "GGG", 0}, //test 2
    {"GGG", "GGGG", 1}, //test 3
    {"GGG", "GGGGG", 0}, //test 4
    {"AAA", "GGG", -6}, //test 5
    {"CCCCCAAGGGGG", "CCCCCGGGGG", 14}, //test 6
    {"TTTTTCCCCCAAGGGGGTTTTT", "TTTTTCCCCCGGGGGTTTTT", 34}, //test 7
    {"AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 40}, //test 8
    {"AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT", -40}, //test 9
    {"AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
     "AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
     160, 2, 2, 5, 1}, //test 10
    {"AAGTGTGTTAATTAATTAATGCTTGTAGGA", "GTTTATGTAGCTTATTCTATCCAAAGCAAT", -12}, //test 11
    {"AAGTGTGTTAATTAATTAATGCTT", "TGTTAATTAATTAATGCTTGGCAAT", 19}, //test 12
    {"GT", "GAT", -1}, //test 13
    {"AAGACATCACGATG", "AAGACACCCCGCACG", 11}, //test 14
    {"GGTT", "GATT", 3, 2, 3, 5, 1}, //test 15
    {"AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAACAAAAAAAAA", 34, 2, 4, 5, 1}, //test 16
    {"GGG", "GGG", 12, 4 /*match*/, 2 /*mismatch*/, 1 /*gap_open*/, 1 /*gap extend*/}, //test 17
    {"GGGGG", "GGGGG", 150, 30, 4, 5, 1}, // test 18
    {"GGGGG", "GGGGG", 250, 50, 2, 5, 1}, // test 19
    {"AAAAA", "AAAA", 2, 2, 4, 6, 1}, // test 20
    {"AAAA", "AAAAA", 2, 2, 4, 6, 1}, // test 21
    {"TTTTT", "TTTT", -1, 0, 1, 1, 1}, // test 22
    {"TTTT", "TTTTT", -1, 0, 1, 1, 1}, // test 23
    {"AAAAAAAAAAAAGAAAAAA",
     "AAAAAAAAAAAAGAAAA",
     27, 2, 4, 6, 1}, // test 24
    {"A", "AAA", -5, 2, 4, 6, 1} // test 25
  };

  // Run all tests if no specific tests are specified
  if (tests_to_run.size() == 0)
  {
    for (int i = 0; i < static_cast<int>(tests.size()); ++i)
      tests_to_run.insert(tests_to_run.end(), i);
  }

  std::string const current_archs = paw::get_current_arch();
  std::cout << "Current archs are: " << current_archs << "\n";
  std::cout << "== global_alignment ==\n";

  auto test_if_expected = [&](long score, Test const & test, long i)
  {
    if (score != test.expected_score)
    {
      std::cout << "\nINCORRECT. Score mismatch in test " << i
                << ". Got score " << score
                << " but I expected " << test.expected_score << "\n\n";
    }
    else
    {
      std::cout << i << " OK!  ";
    }
  };

  for (auto i : tests_to_run)
  {
    assert(i < static_cast<long>(tests.size()));
    auto const & test = tests[i];
    paw::AlignerOptions<uint8_t> opts(false);
    opts.set_match(test.match);
    opts.set_mismatch(test.mismatch);
    opts.set_gap_open(test.gap_open);
    opts.set_gap_extend(test.gap_extend);
    auto score = paw::global_alignment(test.seq1, test.seq2, opts);
    test_if_expected(score, test, i);
  }

  std::cout << "\n== global_alignment_score ==\n";

  for (auto i : tests_to_run)
  {
    auto const & test = tests[i];
    paw::AlignerOptions<uint8_t> opts(false);
    opts.set_match(test.match);
    opts.set_mismatch(test.mismatch);
    opts.set_gap_open(test.gap_open);
    opts.set_gap_extend(test.gap_extend);
    auto score = paw::global_alignment_score(test.seq1, test.seq2, opts);
    test_if_expected(score, test, i);
  }

  std::cout << "\n";
}
