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


inline bool
prefix_matches(std::string const & s1, std::string const & s2)
{
  auto it1 = s1.cbegin();
  auto it2 = s2.cbegin();

  while (it1 != s1.cend() && it2 != s2.cend())
  {
    if (*it1 != *it2)
      return false;

    ++it1;
    ++it2;
  }

  return true;
}


int
main(int argc, char ** argv)
{
  // Local typenames
  //using Ttime = std::chrono::high_resolution_clock;
  //using Tduration = std::chrono::duration<double, std::milli>;

  /// Parse program arguments
  bool backtracking = true;
  std::string fasta_filename;
  std::string fasta_output = "-";
  std::string vcf_output = "-";
  long SPLIT_THRESHOLD = 5;

  try
  {
    paw::Parser parser(argc, argv);
    parser.set_name("Skyr");
    parser.set_version(PAW_VERSION_MAJOR, PAW_VERSION_MINOR, PAW_VERSION_PATCH);
    parser.parse_positional_argument(fasta_filename,
                                     "FASTA",
                                     "A filename of a FASTA file to read sequences from. GZip "
                                     "fasta files supported.");
    parser.parse_option(backtracking,
                        'b',
                        "backtrack",
                        "If set, backtracking/traceback will be made."
                        );
    parser.parse_option(fasta_output,
                        'f',
                        "fasta_output",
                        "Output filename for aligned sequences."
                        );
    parser.parse_option(SPLIT_THRESHOLD, ' ', "split_threshold", "Split threshold to use.");
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

  //auto t0 = Ttime::now();

  // Align the sequence using the SKYR algorithm
  paw::Skyr skyr(fasta.seqs);
  skyr.find_all_edits();
  skyr.find_variants_from_edits();

  //for (auto const & var : skyr.vars)
  //  var.print_seqs(std::cerr);

  auto vars = skyr.split_variants(SPLIT_THRESHOLD);

  for (auto const & var : vars)
    var.print_seqs(std::cerr);

  /*
  auto t1 = Ttime::now();
  std::cout << Tduration(t1 - t0).count() << " ms\n";
  fasta_out.store(fasta_output);

  if (fasta.ids.size() > 0)
  {
    paw::Vcf vcf_out(vcf_output);
    vcf_out.reference = "N" + fasta.seqs[0];

    for (auto const & sn : fasta.ids)
      vcf_out.add_sample_name(sn.substr(1));

    // Shift all positions by 1 for the extra 'N'
    for (auto & var : skyr.vars)
    {
      ++var.pos;

      if (!var.is_snp())
        var.add_base_to_front(vcf_out.reference);

      vcf_out.add_variant(var);
    }

    vcf_out.write();
  }
  */

  return EXIT_SUCCESS;
}
