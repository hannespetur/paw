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
  bool left_column_free = false;
  bool right_column_free = false;

  Test() = default;
  ~Test() = default;

  Test(std::string _seq1,
       std::string _seq2,
       long _expected_score = 0,
       long score_match = 2,
       long score_mismatch = 2,
       long score_gap_open = 5,
       long score_gap_extend = 1,
       bool _left_column_free = false,
       bool _right_column_free = false
       )
    : seq1(_seq1)
    , seq2(_seq2)
    , expected_score(_expected_score)
    , match(score_match)
    , mismatch(score_mismatch)
    , gap_open(score_gap_open)
    , gap_extend(score_gap_extend)
    , left_column_free(_left_column_free)
    , right_column_free(_right_column_free)
  {}
};


template<typename Tuint>
long
calculate_score_from_aligned_strings(paw::AlignmentOptions<Tuint> const & opts,
  std::pair<std::string, std::string> const & a_strings
  )
{
  assert(a_strings.first.size() == a_strings.second.size());
  bool is_ins = false;
  bool is_del = false;
  long score = 0;

  for (long i = 0; i < (long)a_strings.first.size(); ++i)
  {
    auto const & c1 = a_strings.first[i];
    auto const & c2 = a_strings.second[i];

    assert(c1 != '-' || c2 != '-');

    if (c1 == '-')
    {
      // DEL
      if (is_del)
      {
        // DEL extension
        score -= (long)opts.get_gap_extend();
      }
      else
      {
        // DEL open
        score -= (long)opts.get_gap_open();
        is_del = true;
      }

      is_ins = false;
    }
    else if (c2 == '-')
    {
      // INS
      if (is_ins)
      {
        // INS extension
        score -= (long)opts.get_gap_extend();
      }
      else
      {
        // INS open
        score -= (long)opts.get_gap_open();
        is_ins = true;
      }

      is_del = false;
    }
    else
    {
      // Match or mismatch
      if (c1 == c2)
        score += (long)opts.get_match();
      else
        score -= (long)opts.get_mismatch();

      is_del = false;
      is_ins = false;
    }
  }


  if (opts.left_column_free)
  {
    // Add score for every deletion in the beginning
    for (auto it = a_strings.first.begin(); it != a_strings.first.end() && *it == '-'; ++it)
    {
      if (it == a_strings.first.begin())
        score += (long)opts.get_gap_open();
      else
        score += (long)opts.get_gap_extend();
    }
  }

  if (opts.right_column_free)
  {
    // Add score for every deletion in the end
    for (auto it = a_strings.first.rbegin(); it != a_strings.first.rend() && *it == '-'; ++it)
    {
      if (it == a_strings.first.rbegin())
        score += (long)opts.get_gap_open();
      else
        score += (long)opts.get_gap_extend();
    }
  }

  return score;
}


std::vector<std::vector<long> >
transpose(std::vector<std::vector<long> > const & sm)
{
  if (sm.size() == 0)
    return sm;

  std::vector<std::vector<long> > m(sm[0].size(), std::vector<long>(sm.size()));

  for (long i = 0 ; i < (long)sm.size(); ++i)
  {
    for (long j = 0 ; j < (long)sm[i].size(); ++j)
    {
      assert(j < (long)m.size());
      assert(i < (long)m[j].size());
      m[j][i] = sm[i][j];
    }
  }

  return m;
}


/*
void
print_matrix(std::vector<std::vector<long> > const & sm)
{
  for (long i = 0 ; i < (long)sm.size(); ++i)
  {
    for (long j = 0 ; j < (long)sm[i].size(); ++j)
    {
      std::cout << std::setw(6) << sm[i][j];
    }

    std::cout << "\n";
  }
}
*/


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
    {"GGG", "GGG", 6 /*exp. score*/, 2 /*match*/, 2 /*mismatch*/, 10 /*gap_open*/, 1 /*gap extend*/, false /*left column free*/}, //test 0
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
     666, 1, -2, -4, -1}, // test 27
    {"ACGT", "GT", -2, 2, -2, -5, -1}, // test 28
    {"T", "TTTTTCCCCCAAGGGGGTTTTT", -23}, //test 29
    {"GTAGAGGGGGTTGGGCCAAGGTT", "G", -24}, // test 30
    {"GTAGAGGGGGTTGGGCCAAGGTT", "GG", 0, 0, 0, 0, 0}, // test 31
    {"GTAGAGGGGGTTGGGCCAAGGTT", "GTAGGGGGTTGCAGT", 15, 1, 0, 0, 0}, // test 32
    {"GTAGAGGGGGTTGGGCCAAGGTT", "GTAGGGGGTTGCAGT", -8, 0, 1, 1, 1}, // test 33
    {"GGG", "TTTTGGG", 6, 2, -2, -5, -1, true, false}, // test 34
    {"GGTG", "GGTGTCTTGCGTG", 8, 2, -2, -5, -1, false, true}, // test 35
    {"CCCCGTGGGTGGGTGG", "CCCCGGTGGATGGGTGGGGTGTCTTGCGTG", 24, 2, -2, -4, -1, false, true}, // test 36
    {"GGGACGTACGTACGT", "GGCCTTTTGGGACGTACTACGTT", 18, 2, -2, -5, -1, true, false}, // test 37
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

  auto test_if_expected_score = [&](paw::AlignmentOptions<uint8_t> & opts, Test const & test, long i, bool is_swapped) -> bool
  {
    bool are_all_tests_ok = true;
    opts.set_match(test.match).set_mismatch(test.mismatch).set_gap_open(test.gap_open).set_gap_extend(test.gap_extend);
    opts.left_column_free = test.left_column_free;
    opts.right_column_free = test.right_column_free;

    if (is_swapped)
      std::swap(opts.left_column_free, opts.right_column_free);

    paw::AlignmentResults<uint8_t> ar;

    if (is_swapped)
      ar = paw::global_alignment(test.seq2, test.seq1, opts);
    else
      ar = paw::global_alignment(test.seq1, test.seq2, opts);

    std::string test_suffix;

    if (is_swapped)
      test_suffix = " (swapped)";
    else
      test_suffix = " (not swapped)";

    if (ar.score != test.expected_score)
    {
      std::cout << "\nINCORRECT. Final score mismatch in test " << i << test_suffix
                << ". Got score " << ar.score
                << " but I expected " << test.expected_score << "\n" << std::endl;
      are_all_tests_ok = false;
    }

    std::pair<std::string, std::string> aligned_strings;

    if (is_swapped)
      aligned_strings = ar.get_aligned_strings(test.seq2, test.seq1);
    else
      aligned_strings = ar.get_aligned_strings(test.seq1, test.seq2);

    long score_aligned_strings = calculate_score_from_aligned_strings(opts, aligned_strings);

    if (score_aligned_strings != ar.score)
    {
      std::cout << "\nINCORRECT. Traceback error, score from aligned strings does not match in test " << i << test_suffix
                << ". Got score " << score_aligned_strings
                << " but I expected " << test.expected_score << "\n" << std::endl;
      are_all_tests_ok = false;
    }

#ifndef NDEBUG
    assert(opts.score_matrix.size() > 0);
    assert(opts.score_matrix[0].size() > 0);

    if (opts.score_matrix.back().back() != test.expected_score)
    {
      std::cout << "\nINCORRECT. Incorrect final score of the score matrix in test " << i << test_suffix
        << ".  Got score " << opts.score_matrix.back().back()
        << " but I expected " << test.expected_score << "\n" << std::endl;
      are_all_tests_ok = false;
    }
#endif

    // Fix options in case we swapped left_column_free and right_column_free
    if (is_swapped)
      std::swap(opts.left_column_free, opts.right_column_free);

    return are_all_tests_ok;
  };

  for (auto i : tests_to_run)
  {
    assert(i < static_cast<long>(tests.size()));

    auto const & test = tests[i];
    bool is_swapped = !noswap_only && swap_only;
    paw::AlignmentOptions<uint8_t> opts;
    bool are_all_tests_ok = test_if_expected_score(opts, test, i, is_swapped);

    if (!swap_only && !noswap_only && opts.left_column_free == opts.right_column_free)
    {
#ifndef NDEBUG
      std::vector<std::vector<long> > score_matrix = transpose(opts.score_matrix);
      std::vector<std::vector<long> > vE_matrix = transpose(opts.vE_scores);
      std::vector<std::vector<long> > vF_matrix = transpose(opts.vF_scores);

      // Tests for the transpose function
      assert(score_matrix.size() == opts.score_matrix[0].size());
      assert(score_matrix[0].size() == opts.score_matrix.size());
      assert(score_matrix[0][0] == opts.score_matrix[0][0]);
      assert(score_matrix[0][1] == opts.score_matrix[1][0]);
      assert(score_matrix[1][0] == opts.score_matrix[0][1]);
      assert(score_matrix[1][1] == opts.score_matrix[1][1]);
      assert(vE_matrix.size() == opts.vE_scores[0].size());
      assert(vE_matrix[0].size() == opts.vE_scores.size());
      assert(vF_matrix.size() == opts.vF_scores[0].size());
      assert(vF_matrix[0].size() == opts.vF_scores.size());
#endif // NDEBUG

      assert(!is_swapped); // Do align with swapped sequence
      are_all_tests_ok &= test_if_expected_score(opts, test, i, !is_swapped);

#ifndef NDEBUG
      if (score_matrix != opts.score_matrix)
      {
        assert(score_matrix.size() == opts.score_matrix.size());

        for (long r = 0; r < (long) score_matrix.size(); ++r)
        {
          assert(score_matrix[r].size() == opts.score_matrix[r].size());

          for (long c = 0; c < (long) score_matrix[r].size(); ++c)
          {
            if (score_matrix[r][c] != opts.score_matrix[r][c])
            {
              std::cout << "row, col = " << r << "," << c << " mismatch " << score_matrix[r][c] << " != "
                        << opts.score_matrix[r][c] << "\n";
              r = score_matrix.size();
              break;
            }
          }
        }

        are_all_tests_ok = false;
        std::cout << "mismatch in score matrixes\n";
      }

      // Also test if vF and vE are transposed (except for the first row and col)
      assert(vE_matrix.size() == opts.vF_scores.size());

      for (long r = 1; r < static_cast<long>(vE_matrix.size()); ++r)
      {
        assert(vE_matrix[r].size() == opts.vF_scores[r].size());

        for (long c = 1; c < static_cast<long>(vE_matrix[r].size()); ++c)
        {
          if (vE_matrix[r][c] != opts.vF_scores[r][c])
          {
            std::cout << "vE' vs vF. row, col = " << r << "," << c << " mismatch " << vE_matrix[r][c] << " != "
                      << opts.vF_scores[r][c] << "\n";
            r = vE_matrix.size();
            are_all_tests_ok = false;
            break;
          }

          if (vF_matrix[r][c] != opts.vE_scores[r][c])
          {
            std::cout << "vE vs vF'. row, col = " << r << "," << c << " mismatch " << vF_matrix[r][c] << " != "
                      << opts.vE_scores[r][c] << "\n";
            r = vE_matrix.size();
            are_all_tests_ok = false;
            break;
          }
        }
      }
#endif // NDEBUG
    }

    if (are_all_tests_ok)
      ++num_passed_tests;
  }

  std::cout << num_passed_tests << "/" << num_tests << " tests passed." << std::endl;

  if (num_passed_tests != num_tests)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
