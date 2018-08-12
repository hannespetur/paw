#include <paw/align.hpp>
#include <paw/parser.hpp>
#include <paw/internal/config.hpp>

#include <iostream>
#include <string>


int
main(int argc, char ** argv)
{
  std::string seq1;
  std::string seq2;

  try
  {
    paw::Parser parser(argc, argv);
    parser.parse_positional_argument(seq1, "SEQ1", "Sequence 1");
    parser.parse_positional_argument(seq2, "SEQ2", "Sequence 2");
    parser.finalize();
  }
  catch (std::exception const & e)
  {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  paw::AlignmentOptions<uint8_t> opts;
  paw::AlignmentResults<uint8_t> ar = paw::global_alignment(seq1, seq2, opts);
  std::cerr << "Score = " << ar.score << std::endl;
  return EXIT_SUCCESS;
}
