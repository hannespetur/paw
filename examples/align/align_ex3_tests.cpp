#include <paw/align.hpp>
#include <paw/parser.hpp>
#include <paw/internal/config.hpp>

#include <algorithm>
#include <string>
#include <vector>


struct Test
{
  std::string seq1{};
  std::string seq2{};
  long expected_score = 0;
  long match = 2;
  long mismatch = 2;
  long gap_open = 5;
  long gap_extend = 1;

  Test() = default;

  Test(std::string _seq1,
       std::string _seq2,
       long _expected_score = 0,
       long score_match = 2,
       long score_mismatch = 2,
       long score_gap_open = 5,
       long score_gap_extend = 1
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
  bool noswap_only = false;
  bool swap_only = false;

  try
  {
    paw::Parser parser(argc, argv);
    parser.set_name("Alignment example 3 - Tests.");
    parser.parse_option(noswap_only, '1', "noswap", "Set to only run tests with sequences not swapped.");
    parser.parse_option(swap_only, '2', "swap", "Set to run only tests with sequences swapped.");
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

  if (noswap_only && swap_only)
  {
    // Do both
    noswap_only = false;
    swap_only = false;
  }

  std::vector<Test> tests =
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
    {"AAGTGTGTTAATTAATTAATGCTTGTAGGA", "GTTTATGTAGCTTATTCTATCCAAAGCAAT", -12, 2, 2, 5, 1}, //test 11
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
    {"A", "AAA", -5, 2, 4, 6, 1}, // test 25
    {"TGTGTTAATTAATTAATGCTTGTAGGA", "TATGTAGCTTATTCTATCCAAAGCAAT", -6, 2, 2, 5, 1}, //test 26
    {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATA",
     "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA",
     666, 1, -2, -4, -1} // test 27
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
  long num_passed_tests = 0;
  long num_tests = tests_to_run.size();

  if (!swap_only && !noswap_only)
    num_tests *= 2l;

  auto test_if_expected_score = [&](long const score, Test const & test, long i, bool is_swapped)
  {
    std::string test_suffix;

    if (is_swapped)
      test_suffix = " (swap)";

    if (score != test.expected_score)
    {
      std::cout << "\nINCORRECT. Score mismatch in test " << i << test_suffix
                << ". Got score " << score
                << " but I expected " << test.expected_score << "\n" << std::endl;
    }

    ++num_passed_tests;
  };

  for (auto i : tests_to_run)
  {
    assert(i < static_cast<long>(tests.size()));

    auto const & test = tests[i];
    bool is_swapped = !noswap_only && swap_only;
    paw::AlignmentOptions<uint8_t> opts;
    opts.set_match(test.match).set_mismatch(test.mismatch).set_gap_open(test.gap_open).set_gap_extend(test.gap_extend);
    paw::AlignmentResults<uint8_t> ar;

    if (is_swapped)
      ar = paw::global_alignment(test.seq2, test.seq1, opts);
    else
      ar = paw::global_alignment(test.seq1, test.seq2, opts);

    test_if_expected_score(ar.score, test, i, is_swapped);

    if (!swap_only && !noswap_only)
    {
      assert(!is_swapped);
      // Do also swapped
      paw::AlignmentResults<uint8_t> ar_swap = paw::global_alignment(test.seq2, test.seq1, opts);
      test_if_expected_score(ar_swap.score, test, i, !is_swapped);
    }
  }

  std::cout << num_passed_tests << "/" << num_tests << " tests passed." << std::endl;

  //std::cout << "\n== global_alignment_score ==\n";
//
  //for (auto i : tests_to_run)
  //{
  //  auto const & test = tests[i];
  //  paw::AlignmentOptions<uint8_t> opts;
  //  opts.set_match(test.match).set_mismatch(test.mismatch).set_gap_open(test.gap_open).set_gap_extend(test.gap_extend);
  //  paw::AlignmentResults<uint8_t> ar = paw::global_alignment_score(test.seq1, test.seq2, opts);
  //  test_if_expected_score(ar.score, test, i, false);
  //}
//
}
