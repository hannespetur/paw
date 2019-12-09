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
  using Ttime = std::chrono::high_resolution_clock;
  using Tduration = std::chrono::duration<double, std::milli>;

  /// Parse program arguments
  bool backtracking = true;
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

  auto t0 = Ttime::now();

  // Align the sequence using the SKYR algorithm
  paw::SIMDPP_ARCH_NAMESPACE::Skyr skyr(fasta.seqs);
  skyr.find_all_edits();
  skyr.find_variants_from_edits();

  //for (auto const & var : skyr.vars)
  //  var.print_seqs(std::cerr);

  skyr.merge_variants(1);

  /*
  for (std::size_t i = 0; i < skyr.edits.size(); ++i)
  {
    auto const & edit = skyr.edits[i];

    uint32_t del_reach = 0;
    std::string original_seq = skyr.seqs[0]; // Reconstructed sequence
    std::string gapped_seq = skyr.seqs[0]; // Gapped reconstructed sequence
    int64_t shift = 0; // Number of bases shifted by indels
    int64_t gapped_shift = 0; // Number of bases shifted by indels in the gapped sequence

    for (auto & var : skyr.vars)
    {
      assert(var.seqs.size() >= 2);

      auto get_call =
        [del_reach](paw::Variant const & variant, std::set<paw::Event> const & edit)
        {
          std::size_t call = 0;               // Call reference by default
          auto find_it = std::find_if(variant.event_to_allele.begin(),
                                      variant.event_to_allele.end(),
                                      [&](std::pair<paw::Event, uint32_t> const & e) -> bool
        {
          return edit.count(e.first);
        });

          // Check if the sequence had any event at this variant location
          if (find_it != variant.event_to_allele.end())
          {
            call = find_it->second;
          }
          else if (variant.pos < del_reach)
          {
            // Call asterisk if the variant is deleted by a previous deletion
            assert(variant.seqs.back() == "*");
            call = variant.seqs.size() - 1;
          }

          return call;
        };

      uint16_t const call = get_call(var, edit);
      var.add_call(call);

      // Check if we need to replace strings in the reconstructed sequence
      if (call != 0 && var.seqs[call] != "*")
      {
        int64_t new_pos = var.pos + shift;
        original_seq = original_seq.replace(new_pos, var.seqs[0].size(), var.seqs[call]);
        shift += var.seqs[call].size() - var.seqs[0].size();
      }

      if (var.seqs[call] != "*")
      {
        int64_t new_gapped_pos = var.pos + gapped_shift;
        auto longest_seq = var.seqs.cbegin();

        if (call != 0)
          gapped_seq = gapped_seq.replace(new_gapped_pos, var.seqs[0].size(), var.seqs[call]);

        if (var.is_insertion())
        {
          // The longest sequence is at the back
          longest_seq = var.seqs.cbegin() + (var.seqs.size() - 1);

          // Add gaps at end
          if (longest_seq->size() > var.seqs[call].size())
          {
            std::string gaps(longest_seq->size() - var.seqs[call].size(), '-');
            gapped_seq = gapped_seq.replace(new_gapped_pos + var.seqs[call].size(), 0, gaps);
          }

          gapped_shift += longest_seq->size() - var.seqs[0].size();
        }
        else if (var.is_deletion())
        {
          // Longest sequence is the reference
          if (longest_seq->size() > var.seqs[call].size())
          {
            std::string gaps(longest_seq->size() - var.seqs[call].size(), '-');
            gapped_seq = gapped_seq.replace(new_gapped_pos, 0, gaps);
          }
        }
      }

      if (var.is_deletion())
      {
        assert(var.seqs[0].size() >= var.seqs[call].size());
        del_reach = std::max(del_reach,
                             static_cast<uint32_t>(var.pos +
                                                   var.seqs[0].size() -
                                                   var.seqs[call].size()
                                                   )
                             );
      }
    }

    fasta_out.add_record(fasta.ids[i], gapped_seq);
  }
  */

  auto t1 = Ttime::now();
  std::cout << Tduration(t1 - t0).count() << " ms\n";
  fasta_out.store(fasta_output);

  if (fasta.ids.size() > 0)
  {
    paw::SIMDPP_ARCH_NAMESPACE::Vcf vcf_out(vcf_output);
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
