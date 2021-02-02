#include <paw/align.hpp>
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
  opts.get_aligned_strings = true;

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
    auto const & id = contigs.ids[s];
    auto const & seq = contigs.seqs[s];

    if (static_cast<long>(seq.size()) < min_length)
      continue;

    paw::global_alignment(seq, ref.seqs[0], opts);
    auto ar1 = opts.get_alignment_results();
    auto aligned_strings1 = *(ar1->aligned_strings_ptr);

    std::string rseq = paw::reverse_complement(seq);
    paw::global_alignment(rseq, ref.seqs[0], opts);
    auto ar2 = opts.get_alignment_results();
    auto aligned_strings2 = *(ar2->aligned_strings_ptr);

    if (ar2->score > ar1->score)
    {
      // make aligned_strings1 better
      std::swap(aligned_strings1, aligned_strings2);
    }

    std::swap(aligned_strings1.first, aligned_strings1.second);
    auto edits = paw::get_edit_script(aligned_strings1, true, false); // is_normalize, is_trim_indel_on_ends

    if (edits.size() > 0)
    {
      //std::cerr << "Found " << edits.size() << " variants.\n";
      bool is_id_printed = false;
      //std::cout << id << "\n";
      std::vector<paw::Event2> edits_vec(edits.begin(), edits.end());
      //auto it = edits.begin();

      for (long i = 0; i < static_cast<long>(edits.size()); ++i)
      {
        //assert(it != edits.end());
        auto const & e = edits_vec[i];

        if (!e.is_snp() &&
            (i == 0 || i == static_cast<long>(edits_vec.size() - 1) || static_cast<long>(e.pos + e.ref.size()) > 29850))
        {
          continue;
        }

        long const size_diff = static_cast<long>(e.ref.size()) - static_cast<long>(e.alt.size());
        //long const threshold = 30;

        // Don't allow deletions to be larger than 1k bp, since that it larger than primers
        if (size_diff >= 500)
          continue;

        if ((i >= 2 && (e.pos - 100) < edits_vec[i - 2].pos) ||
            (i <= static_cast<long>(edits_vec.size() - 3) && (e.pos + 100) > edits_vec[i + 2].pos))
        {
          continue;
        }

        //if (size_diff < -threshold || size_diff > threshold)
        {
          if (!is_id_printed)
          {
            std::cout << id << "\n";
            is_id_printed = true;
          }

          std::cout << (e.pos + 1) << " " << (e.ref.size() == 0 ? "-" : e.ref) << " "
                    << (e.alt.size() == 0 ? "-" : e.alt) << " "
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
