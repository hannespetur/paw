#include <paw/align.hpp>
#include <paw/internal/config.hpp>

#include <chrono> // std::chrono::high_resolution_clock
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
main(int, char **)
{
  std::vector<Test> const tests = {
      {"GGG", "GGG", 6},
      {"GGGG", "GGG", 1},
      {"GGGGG", "GGG", 0},
      {"GGG", "GGGG", 1},
      {"GGG", "GGGGG", 0},
      {"AAA", "GGG", -6},
      {"CCCCCAAGGGGG", "CCCCCGGGGG", 14},
      {"TTTTTCCCCCAAGGGGGTTTTT", "TTTTTCCCCCGGGGGTTTTT", 34},
      {"AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 40},
      {"AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT", -40},
      {"AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
       "AAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAA",
       160
      },
      {"AAGTGTGTTAATTAATTAATGCTTGTAGGA", "GTTTATGTAGCTTATTCTATCCAAAGCAAT", -12},
      {"AAGTGTGTTAATTAATTAATGCTT", "TGTTAATTAATTAATGCTTGGCAAT", 19},
      {"GT", "GAT", -1},
      {"AAGACATCACGATG", "AAGACACCCCGCACG", 11}
  };

  std::string const current_archs = paw::get_current_arch();
  std::cerr << "Current archs are: " << current_archs << "\n";

  for (long i = 0; i < static_cast<long>(tests.size()); ++i)
  {
    auto const & test = tests[i];
    paw::AlignerOptions<uint8_t> opts;
    auto score = paw::align(test.seq1, test.seq2, opts);

    if (score != test.expected_score)
    {
      std::cerr << "Score mismatch in test " << i
                << ". Got score " << score
                << " but I expected " << test.expected_score << std::endl;
    }
  }
}
