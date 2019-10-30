#include <chrono>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <set>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <paw/parser.hpp>
#include <paw/internal/config.hpp>

#include <paw/align.hpp>


namespace io = boost::iostreams;


int
main(int argc, char ** argv)
{
  std::string fasta_filename;
  std::string fasta_output = "-";
  std::string vcf_output = "-";

  try
  {
    paw::Parser parser(argc, argv);
    parser.set_name("Skyr");
    parser.set_version(PAW_VERSION_MAJOR, PAW_VERSION_MINOR, PAW_VERSION_PATCH);
    parser.parse_positional_argument(fasta_filename,
                                     "FASTA",
                                     "A filename of a FASTA file to read sequences from. GZip "
                                     "fasta files supported.");
    parser.parse_option(fasta_output,
                        'f',
                        "fasta_output",
                        "Output filename for aligned sequences."
                        );
    parser.parse_option(vcf_output, 'o', "vcf_output", "Output filename with variant records.");
    parser.finalize();
  }
  catch (const std::exception & e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }
  ///

  // Load records from the FASTA file
  paw::Fasta fasta;
  fasta.load(fasta_filename);
  paw::Fasta fasta_out;

  // If less then two records, there is nothing to do
  if (fasta.seqs.size() < 2)
  {
    std::cerr << "[skyr] WARNING: Cannot align, only " << fasta.seqs.size()
              << " sequences found in " << fasta_filename << "\n";
    return EXIT_SUCCESS;
  }


  /// print sequences
  {
    std::cout << "Sequences:\n";

    for (auto const & seq : fasta.seqs)
    {
      std::cout << seq << "\n";
    }
  }
  ///

  paw::SIMDPP_ARCH_NAMESPACE::Skyr skyr(fasta.seqs);
  skyr.find_all_edits();
  std::cout << "Total number of edits: " << skyr.all_edits.size() << "\n";
  skyr.find_variants_from_edits();

  /// print variants founds from the edits
  {
    std::cout << "Num edits: " << skyr.edits.size() << "\n";
    std::cout << "Num variants: " << skyr.vars.size() << "\n";

    for (auto const & variant : skyr.vars)
    {
      variant.print_seqs(std::cout);
    }
  }

  return EXIT_SUCCESS;
}
