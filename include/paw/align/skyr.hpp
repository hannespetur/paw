#pragma once

#include <set> // std::set<T>
#include <string> // std::string
#include <vector> // std::vector<T>

#include <paw/align/global_alignment.hpp>
#include <paw/align/event.hpp>
#include <paw/align/variant.hpp>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

class Skyr
{
public:
  using Tscores = std::vector<int64_t>;   // Container for alignment scores

  std::vector<std::string> seqs;
  std::set<Event2> free_edits;   // Alignment edits that have been made free
  std::vector<std::set<Event2> > edits;   // Edits found for each alignment
  std::multiset<Event2> all_edits;
  std::vector<Variant> vars;

  Skyr(std::vector<std::string> const & _seqs);
  Skyr(std::vector<std::vector<char> > const & _seqs);

  void find_all_edits();
  void find_variants_from_edits();

private:
  std::size_t find_most_similar_haplotype(Tscores const & scores) const;
};


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw


#if defined(IMPLEMENT_PAW) || defined(__JETBRAINS_IDE__)


#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

//#include "boost_simd_align.hpp"


namespace
{

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


} // anon namespace


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{

Skyr::Skyr(std::vector<std::string> const & _seqs)
  : seqs(_seqs)
  , edits(_seqs.size())
{}

Skyr::Skyr(std::vector<std::vector<char> > const & _seqs)
  : edits(_seqs.size())
{
  seqs.reserve(_seqs.size());

  for (auto const & _seq : _seqs)
    seqs.push_back(std::string(_seq.cbegin(), _seq.cend()));
}


std::size_t
Skyr::find_most_similar_haplotype(Tscores const & scores) const
{
  int64_t max_score = std::numeric_limits<int64_t>::min();
  std::size_t max_events = static_cast<std::size_t>(-1);
  std::size_t max_events_seen = 0;
  std::size_t max_i = static_cast<std::size_t>(-1);

  for (std::size_t i = 1; i < edits.size(); ++i)
  {
    std::size_t event_seen_count = 0;

    auto is_novel = [&](Event2 const & e)
                    {
                      return e.is_snp() && free_edits.count(e) == 0;
                    };

    if (std::none_of(edits[i].begin(), edits[i].end(), is_novel))
      continue;

    for (auto const & e : edits[i])
      event_seen_count += all_edits.count(e);

    if (scores[i] > max_score ||
        (scores[i] == max_score && edits[i].size() < max_events) ||
        (scores[i] == max_score && edits[i].size() == max_events &&
         event_seen_count > max_events_seen)
        )
    {
      max_score = scores[i];
      max_events = edits[i].size();
      max_events_seen = event_seen_count;
      max_i = i;
    }
  }

  return max_i;
}


void
Skyr::find_all_edits()
{
  std::size_t iteration = 0;
  Tscores scores(seqs.size());

  while (true)
  {
    ++iteration;
    std::cout << "Iteration #" << iteration << "\n";
    using Tuint = uint8_t;
    AlignmentOptions<Tuint> opts;
    //opts.backtracking = true; // Force backtracking

    // Construct an aligner with the reference sequence
    //paw::global_alignment(seqs[0].cbegin(), seqs[0].cend(), opts); //, free_edits);
    all_edits.clear();

    for (std::size_t i = 1; i < seqs.size(); ++i)
    {
      global_alignment<std::string, Tuint>(seqs[0], seqs[i], opts);
      auto ar = opts.get_alignment_results();
      //auto ac = opts.get_alignment_cache();
      scores[i] = ar->score;
      std::cerr << "Score = " << ar->score << "\n";

      //scores[i] = aligner.align(seqs[i].cbegin(), seqs[i].cend());
      auto aligned_strings = ar->get_aligned_strings(seqs[0], seqs[i]);
      std::cout << aligned_strings.first << "\n"
                << aligned_strings.second << "\n";
      //edits[i] = aligner.get_edit_script(aligned_strings);
      //all_edits.insert(edits[i].begin(), edits[i].end());
    }

    std::size_t const max_i = find_most_similar_haplotype(scores);

    if (max_i != static_cast<std::size_t>(-1))
    {
      for (auto const & e : edits[max_i])
      {
        if (e.is_snp())
          free_edits.insert(e);
      }
    }
    else
    {
      // We are finished, no new edits
      break;
    }
  }
}


void
Skyr::find_variants_from_edits()
{
  Variant new_var;

  for (auto it = all_edits.begin(); it != all_edits.end(); ++it)
  {
    {
      auto next_it = std::next(it);

      if (next_it != all_edits.end() && *it == *next_it)
        continue;   // Same as next edit, we can skip it
    }

    auto const & e = *it;

    if (!new_var.has_sequences())
    {
      new_var = {static_cast<uint32_t>(e.pos), {e.ref}};
      new_var.add_event(e);
      continue;
    }

    if (e.pos == new_var.pos && e.is_insertion() && new_var.is_insertion() &&
        prefix_matches(new_var.seqs.back(), e.alt)
        )
    {
      // Case 1: insertions at the same position
      new_var.add_event(e);
    }
    else if (e.pos == new_var.pos && e.is_deletion() && new_var.is_deletion())
    {
      // Case 2: Deletions at the same position
      assert(e.ref.size() > new_var.seqs[0].size());

      /// Extend all sequences
      auto begin_extend = e.ref.cbegin() + new_var.seqs[0].size();

      for (auto & s : new_var.seqs)
        s.append(begin_extend, e.ref.cend());
      ///

      new_var.add_event(e);
    }
    else if (e.pos == new_var.pos && e.ref.size() == 1 && e.alt.size() == 1 &&
             new_var.seqs[0].size() == 1 && new_var.seqs[1].size() == 1
             )
    {
      // Case 3: both are SNPs at the same position
      new_var.add_event(e);
    }
    else
    {
      // Otherwise, they are two different events
      vars.push_back(std::move(new_var));
      new_var = {static_cast<uint32_t>(e.pos), {e.ref} };
      new_var.add_event(e);
    }
  }

  if (new_var.seqs.size() > 1)
    vars.push_back(std::move(new_var));

  // Add asterisk where needed
  {
    uint32_t del_start = 0;
    uint32_t del_reach = 0;

    for (auto & v : vars)
    {
      if (v.pos > del_start && v.pos + v.seqs[0].size() <= del_reach)
        v.seqs.push_back("*");

      del_start = v.pos;
      del_reach = std::max(del_reach, static_cast<uint32_t>(v.pos + v.seqs[0].size()));
    }
  }
}


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw

#endif // IMPLEMENT_PAW
