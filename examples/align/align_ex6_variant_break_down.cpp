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

  // If less then two records, there is nothing to do
  if (fasta.seqs.size() < 2)
  {
    std::cerr << "[skyr] WARNING: Cannot align, only " << fasta.seqs.size()
              << " sequences found in " << fasta_filename << "\n";
    return EXIT_SUCCESS;
  }

  paw::Skyr skyr(fasta.seqs);
  skyr.find_all_edits(true); // normalize
  skyr.find_variants_from_edits();
  skyr.populate_variants_with_calls();

  for (auto & var : skyr.vars)
  {
    std::cout << var.pos << '\t';

    // Loop over alleles
    for (long s = 0; s < static_cast<long>(var.seqs.size()); ++s)
    {
      if (s > 0)
        std::cout << ",";

      if (var.seqs[s].size() == 0)
        std::cout << '-';
      else
        std::cout << var.seqs[s];
    }

    std::cout << '\t';

    // Loop over alignments
    for (long i = 0; i < static_cast<long>(skyr.seqs.size()); ++i)
    {
      if (i > 0)
        std::cout << '\t';

      std::cout << std::to_string(var.get_call(i));
    }

    std::cout << '\n';
  }

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

  return EXIT_SUCCESS;
}
