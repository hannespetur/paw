//#include <paw/align/event.hpp>
#include <paw/align.hpp>
#include <paw/align/sequence_utils.hpp>
#include <paw/parser.hpp>
#include <paw/internal/config.hpp>

#include <iostream>
#include <string>


int
main(int argc, char ** argv)
{
  std::string ref_fn;
  std::string contigs_fn;
  int min_length = 300;

  try
  {
    paw::Parser parser(argc, argv);
    parser.parse_positional_argument(ref_fn, "REF", "Reference FASTA sequence.");
    parser.parse_positional_argument(contigs_fn, "FASTA", "Contig FASTA with sequences to align to the ref.");
    parser.parse_option(min_length, 'l', "min_length", "Minimum length for a contig for it to be considered.");
    parser.finalize();
  }
  catch (std::exception const & e)
  {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  paw::AlignmentOptions<uint8_t> opts;
  opts.left_column_free = true;
  opts.right_column_free = true;

  paw::Fasta ref;
  ref.load(ref_fn);

  if (ref.seqs.size() != 1)
  {
    std::cerr << "Expected only one reference contig, got " << ref.seqs.size() << "\n";
    std::exit(1);
  }

  paw::Fasta contigs;
  contigs.load(contigs_fn);
  //std::cerr << "Ref num seqs = " << contigs.seqs.size() << "\n";

  for (long s = 0; s < static_cast<long>(contigs.seqs.size()); ++s)
  {
    auto const & seq = contigs.seqs[s];

    if (static_cast<long>(seq.size()) < min_length)
      continue;

    paw::global_alignment(seq, ref.seqs[0], opts);
    auto ar1 = opts.get_alignment_results();
    auto aligned_strings1 = ar1->get_aligned_strings(seq, ref.seqs[0]);

    std::string rseq = paw::reverse_complement(seq);
    paw::global_alignment(rseq, ref.seqs[0], opts);
    auto ar2 = opts.get_alignment_results();
    auto aligned_strings2 = ar1->get_aligned_strings(rseq, ref.seqs[0]);

    if (ar2->score > ar1->score)
    {
      // make aligned_strings1 better
      std::swap(aligned_strings1, aligned_strings2);
    }

    std::swap(aligned_strings1.first, aligned_strings1.second);
    //std::cerr << aligned_strings.first << "\n" << aligned_strings.second << "\n";
    auto edits = paw::get_edit_script(aligned_strings1, true); // is_normalize

    if (edits.size() > 0)
    {
      //std::cerr << aligned_strings.first << "\n";
      //std::cerr << "Found " << edits.size() << " variants.\n";
      auto it = edits.begin();

      for (long i = 0; i < static_cast<long>(edits.size()); ++i, ++it)
      {
        assert(it != edits.end());
        auto const & e = *it;

        if (e.is_deletion() &&
            (i == 0 || i == static_cast<long>(edits.size() - 1) || static_cast<long>(e.pos + e.ref.size()) > 29850))
        {
          continue;
        }

        long const size_diff = static_cast<long>(e.ref.size()) - static_cast<long>(e.alt.size());
        long const threshold = 30;

        // Don't allow deletions to be larger than 1k bp, since that it larger than primers
        if (size_diff >= 1000)
          continue;

        if (size_diff < -threshold || size_diff > threshold)
        {
          std::cout << (e.pos + 1) << " " << e.ref << " " << e.alt << " "
                    << e.ref.size() << " " << e.alt.size() << "\n";
        }
      }
    }

    //std::cerr << "Score = " << ar->score << std::endl;
  }
  //paw::global_alignment(seq1, seq2, opts);
  //std::cerr << "Score = " << opts.get_alignment_results()->score << std::endl;
  return EXIT_SUCCESS;
}
