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
  long match{1};
  long mismatch{4};
  long gap_open{6};
  long gap_extend{1};

  Test() = default;
  ~Test() = default;

  Test(std::string const & _seq1,
       std::string const & _seq2,
       long _expected_score = 0,
       long score_match = 1,
       long score_mismatch = 4,
       long score_gap_open = 6,
       long score_gap_extend = 1)
    : seq1(_seq1)
    , seq2(_seq2)
    , expected_score(_expected_score)
    , match(score_match)
    , mismatch(score_mismatch)
    , gap_open(score_gap_open)
    , gap_extend(score_gap_extend)
  {}
};


template <typename Tuint>
long
calculate_score_from_aligned_strings(paw::AlignmentOptions<Tuint> const & opts,
                                     std::pair<std::string, std::string> const & a_strings
                                     )
{
  assert(a_strings.first.size() == a_strings.second.size());
  bool is_ins{false};
  bool is_del{false};
  long score{0};

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
      if (c1 == c2 || c1 == 'N' || c2 == 'N')
        score += (long)opts.get_match();
      else
        score -= (long)opts.get_mismatch();

      is_del = false;
      is_ins = false;
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

  for (long i = 0; i < (long)sm.size(); ++i)
  {
    for (long j = 0; j < (long)sm[i].size(); ++j)
    {
      assert(j < (long)m.size());
      assert(i < (long)m[j].size());
      m[j][i] = sm[i][j];
    }
  }

  return m;
}


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


int
main(int argc, char ** argv)
{
  using Tuint = uint8_t;
  std::vector<int> tests_to_run;

  try
  {
    paw::Parser parser(argc, argv);
    parser.set_name("Alignment example 3 - Tests.");
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

  std::vector<Test> tests =
  {
    {"G", "G", 1, 1, 4, 6, 1}, // test 0
    {"GG", "GG", 2, 1, 4, 6, 1}, // test 1
    {"GG", "GA", -3, 1, 4, 6, 1}, // test 2
    {"GGG", "GAA", -4, 1, 4, 6, 1}, // test 3
    {"GGAGG", "GGGGG", 0, 1, 4, 6, 1}, // test 4
    {"GGAAA", "GGGGG", -3, 1, 4, 6, 1}, // test 5
    {"GGAAA", "GGGGGGGGGGGGGGG", -3, 1, 4, 6, 1}, // test 6
    {"GGAAA", "GGGGGGGGGGGGGGGA", -3, 1, 4, 6, 1}, // test 7
    {"GGAAA", "GGAAA", 5, 1, 4, 6, 1}, // test 8
    {"GGAAA", "GGGAAA", 3, 1, 2, 2, 1}, // test 9
    {"GGAAA", "GGGAAA", 3, 1, 4, 2, 1}, // test 10
    {"GGAAA", "GGGAAA", 2, 1, 4, 3, 1}, // test 11
    {"GGAAA", "GGGAAA", 2, 1, 5, 3, 1}, // test 12
    {"GGAAA", "GGGAAA", 2, 1, 6, 3, 1}, // test 13
    {"GGAAA", "GGGAAA", 2, 1, 7, 3, 1}, // test 14
    {"GGAAA", "GGGAAA", 1, 1, 3, 4, 1}, // test 15
    {"GGAAA", "GGGAAA", 0, 1, 4, 5, 1}, // test 16
    {"GGAAA", "GGGAAA", 1, 1, 3, 5, 1}, // test 17
    {"GGAAA", "GGGAAAAAAAAAAAAA", 1, 1, 3, 5, 1}, // test 18
    {"GGAAA", "GGGAAAAAAAAAAAAA", 0, 1, 4, 6, 1}, // test 19
    {"GGCCCCCCCCCCCCC", "GGGAAAAAAAAAAAAAAAAA", -3, 1, 4, 6, 1}, // test 19
    {"GGCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "GGCCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", 25, 1, 4, 6, 1}, // test 20
    {"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", 60, 1, 4, 6, 1}, // test 21
    {"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", 60, 1, 8, 12, 1}, // test 22
    {"GGAAA", "GGGAAAAAAAAAAAAA", 4, 2, 4, 6, 1}, // test 23
    {"GGAAA", "GGGAAAAAAAAAAAAA", 2, 1, 2, 6, 1}, // test 24
    {"GGGAAA", "GGAAA", 0, 1, 4, 5, 1}, // test 25
    {"GGGAAA", "GGAAA", 0, 1, 3, 5, 1}, // test 26
    {"GGGAAA", "GGAAACCCCCCC", -1, 1, 4, 6, 1}, // test 27
    {"GGAAA", "GGGAAACCCCCCC", 0, 1, 4, 6, 1}, // test 28

//     {"GGGG", "GGG", 1}, // test 1
//     {"GGGGG", "GGG", 3}, // test 2
//     {"GGG", "GGGG", 1}, // test 3
//     {"GGG", "GGGGG", 0}, // test 4
//     {"AAA", "GGG", -6}, // test 5
//     {"CCCCCAAGGGGG", "CCCCCGGGGG", 14}, // test 6
//     {"TTTTTCCCCCAAGGGGGTTTTT", "TTTTTCCCCCGGGGGTTTTT", 34}, // test 7
//     {"AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 40}, //test 8
//     {"AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT", -40}, //test 9
//     {"AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
//      "AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
//      160, 2, 2, 5, 1}, //test 10
//     {"AAGTGTGTTAATTAATTAATGCTTGTAGGA", "GTTTATGTAGCTTATTCTATCCAAAGCAAT", -12, 2, 2, 5, 1}, //test 11
//     {"AAGTGTGTTAATTAATTAATGCTT", "TGTTAATTAATTAATGCTTGGCAAT", 19}, //test 12
//     {"GT", "GAT", -1}, //test 13
//     {"AAGACATCACGATG", "AAGACACCCCGCACG", 11}, //test 14
//     {"GGTT", "GATT", 3, 2, 3, 5, 1}, //test 15
//     {"AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAACAAAAAAAAA", 34, 2, 4, 5, 1}, //test 16
//     {"GGG", "GGG", 12, 4 /*match*/, 2 /*mismatch*/, 1 /*gap_open*/, 1 /*gap extend*/}, //test 17
//     {"GGGGG", "GGGGG", 150, 30, 4, 5, 1}, // test 18
//     {"GGGGG", "GGGGG", 200, 40, 2, 5, 1}, // test 19
//     {"AAAAA", "AAAA", 2, 2, 4, 6, 1}, // test 20
//     {"AAAA", "AAAAA", 2, 2, 4, 6, 1}, // test 21
//     {"TTTTT", "TTTT", -1, 0, 1, 1, 1}, // test 22
//     {"TTTT", "TTTTT", -1, 0, 1, 1, 1}, // test 23
//     {"AAAAAAAAAAAAGAAAAAA",
//      "AAAAAAAAAAAAGAAAA",
//      27, 2, 4, 6, 1}, // test 24
//     {"A", "AAA", -5, 2, 4, 6, 1}, // test 25
//     {"TGTGTTAATTAATTAATGCTTGTAGGA", "TATGTAGCTTATTCTATCCAAAGCAAT", -6, 2, 2, 5, 1}, //test 26
//     {
//       "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATATCTATATATATATACATATATATATA",
//       "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA",
//       666, 1, -2, -4, -1}, // test 27
//     {"ACGT", "GT", -2, 2, -2, -5, -1}, // test 28
//     {"T", "TTTTTCCCCCAAGGGGGTTTTT", -23}, //test 29
//     {"GTAGAGGGGGTTGGGCCAAGGTT", "G", -24}, // test 30
//     {"GTAGAGGGGGTTGGGCCAAGGTT", "GG", 0, 0, 0, 0, 0}, // test 31
//     {"GTAGAGGGGGTTGGGCCAAGGTT", "GTAGGGGGTTGCAGT", 15, 1, 0, 0, 0}, // test 32
//     {"GTAGAGGGGGTTGGGCCAAGGTT", "GTAGGGGGTTGCAGT", -8, 0, 1, 1, 1}, // test 33
  };

  // Run all tests if no specific tests are specified
  if (tests_to_run.size() == 0)
  {
    for (int i = 0; i < static_cast<int>(tests.size()); ++i)
      tests_to_run.insert(tests_to_run.end(), i);
  }

  std::string const current_archs = paw::get_current_arch();
  std::cout << "Current archs are: " << current_archs << "\n";
  std::cout << "== pairwise_ext_alignment ==\n";
  long num_passed_tests{0};
  long num_tests = tests_to_run.size();

  auto test_if_expected_score =
    [&](paw::AlignmentOptions<Tuint> & opts, Test const & test, long i) -> bool
    {
      bool are_all_tests_ok = true;
      opts.set_match(test.match).set_mismatch(test.mismatch);
      opts.set_gap_open(test.gap_open).set_gap_extend(test.gap_extend);
      opts.get_aligned_strings = true;

      paw::pairwise_ext_alignment(test.seq1, test.seq2, opts);

      paw::AlignmentResults const & ar = *opts.get_alignment_results();
      std::string test_suffix;

      if (ar.score != test.expected_score)
      {
        std::cout << "\nINCORRECT. Final score mismatch in test " << i << test_suffix
                  << ". Got score " << ar.score
                  << " but I expected " << test.expected_score << "\n" << std::endl;
        are_all_tests_ok = false;
      }

      std::pair<std::string, std::string> const & aligned_strings = *ar.aligned_strings_ptr;
      long score_aligned_strings = calculate_score_from_aligned_strings(opts, aligned_strings);

      if (ar.query_end != static_cast<int>(test.seq1.size()))
        score_aligned_strings -= opts.get_clip();

      if (score_aligned_strings != test.expected_score)
      {
        std::cout << "\nINCORRECT. Traceback error, score from aligned strings does not match in test "
                  << i << test_suffix
                  << ". Got score " << score_aligned_strings
                  << " but I expected " << test.expected_score << "\n" << std::endl;

        std::cout << aligned_strings.first << '\n' << aligned_strings.second << '\n';

        are_all_tests_ok = false;
      }

      /*
#ifndef NDEBUG
      assert(opts.score_matrix.size() > 0);
      assert(opts.score_matrix[0].size() > 0);

      if (opts.score_matrix.back().back() != test.expected_score)
      {
        std::cout << "\nINCORRECT. Incorrect final score of the score matrix in test "
                  << i << test_suffix
                  << ".  Got score " << opts.score_matrix.back().back()
                  << " but I expected " << test.expected_score << "\n" << std::endl;

        are_all_tests_ok = false;
      }
#endif
*/

      return are_all_tests_ok;
    };

  for (auto i : tests_to_run)
  {
    assert(i < static_cast<long>(tests.size()));

    auto const & test = tests[i];
    paw::AlignmentOptions<Tuint> opts;
    bool are_all_tests_ok = test_if_expected_score(opts, test, i);

      are_all_tests_ok &= test_if_expected_score(opts, test, i);

      /*
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
      assert(vF_matrix.size() == opts.vE_scores.size());

      for (long r = 1; r < static_cast<long>(vE_matrix.size()); ++r)
      {
        assert(vE_matrix[r].size() == opts.vF_scores[r].size());
        assert(vF_matrix[r].size() == opts.vE_scores[r].size());
        assert(vF_matrix[r].size() == vE_matrix[r].size());

        for (long c = 1; c < static_cast<long>(vE_matrix[r].size()); ++c)
        {
          if (vE_matrix[r][c] != opts.vF_scores[r][c])
          {
            std::cout << "INCORRECT: Test " << i << ": vE' vs vF. row, col = " << r << "," << c << " mismatch " <<
              vE_matrix[r][c] << " != "
                      << opts.vF_scores[r][c] << "\n";
            r = vE_matrix.size();
            are_all_tests_ok = false;
            break;
          }

          if (vF_matrix[r][c] != opts.vE_scores[r][c])
          {
            std::cout << "INCORRECT: Test " << i << " vE vs vF'. row, col = " << r << "," << c << " mismatch " <<
              vF_matrix[r][c] << " != "
                      << opts.vE_scores[r][c] << "\n";
            r = vE_matrix.size();
            are_all_tests_ok = false;
            break;
          }
        }
      }
#endif // NDEBUG
      */

    if (are_all_tests_ok)
      ++num_passed_tests;
  }

  std::cout << num_passed_tests << "/" << num_tests << " tests passed." << std::endl;

  if (num_passed_tests != num_tests)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
