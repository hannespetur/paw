#pragma once

#include <cstdint> // uint8_t
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
  std::vector<uint8_t> is_done;

  Skyr(std::vector<std::string> const & _seqs);
  Skyr(std::vector<std::vector<char> > const & _seqs);

  void find_all_edits();
  void find_variants_from_edits();
  void populate_variants_with_calls();

private:
  long find_most_similar_haplotype(Tscores const & scores) const;
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
  , is_done(_seqs.size(), 0)
{}


Skyr::Skyr(std::vector<std::vector<char> > const & _seqs)
  : edits(_seqs.size())
  , is_done(_seqs.size(), 0)
{
  seqs.reserve(_seqs.size());

  for (auto const & _seq : _seqs)
    seqs.push_back(std::string(_seq.cbegin(), _seq.cend()));
}


long
Skyr::find_most_similar_haplotype(Tscores const & scores) const
{
  long max_score = std::numeric_limits<long>::min();
  long max_events = -1;
  long max_events_seen = 0;
  long max_i = -1;

  for (long i = 1; i < static_cast<long>(edits.size()); ++i)
  {
    assert(i < static_cast<long>(is_done.size()));

    if (is_done[i] != 0)
      continue;

    long event_seen_count = 0;

    for (auto const & e : edits[i])
      event_seen_count += all_edits.count(e);

    if (scores[i] > max_score ||
        (scores[i] == max_score && static_cast<long>(edits[i].size()) < max_events) ||
        (scores[i] == max_score && static_cast<long>(edits[i].size()) == max_events &&
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
  Tscores scores(seqs.size(), std::numeric_limits<long>::min());
  using Tuint = uint8_t;
  AlignmentOptions<Tuint> opts;

  while (std::find(is_done.begin() + 1, is_done.end(), 0) != is_done.end())
  {
    ++iteration;
    all_edits.clear();

    for (long i = 1; i < static_cast<long>(seqs.size()); ++i)
    {
      assert(i < static_cast<long>(is_done.size()));
      assert(i < static_cast<long>(edits.size()));

      if (is_done[i] == 1)
      {
        all_edits.insert(edits[i].begin(), edits[i].end());
        continue;
      }

      global_alignment<std::string, Tuint>(seqs[0], seqs[i], opts);
      auto ar = opts.get_alignment_results();
      scores[i] = ar->score;

      auto aligned_strings = ar->get_aligned_strings(seqs[0], seqs[i]);
      edits[i] = get_edit_script(aligned_strings);
      all_edits.insert(edits[i].begin(), edits[i].end());
    }

    long max_i = find_most_similar_haplotype(scores);

    while (max_i != -1)
    {
      assert(max_i < static_cast<long>(seqs.size()));
      is_done[max_i] = 1;
      bool is_novel_snp = false;

      for (auto const & e : edits[max_i])
      {
        if (e.is_snp() && free_edits.count(e) == 0)
        {
          is_novel_snp = true;
          free_edits.insert(e);
          opts.get_alignment_cache()->set_free_snp(e.pos, e.alt[0]);
        }
      }

      if (is_novel_snp)
        break; // We added a novel SNP, time to remap.
      else
        max_i = find_most_similar_haplotype(scores);
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

    if (e.pos == new_var.pos && e.is_insertion() && new_var.is_insertion())
// && prefix_matches(new_var.seqs.back(), e.alt)
    {
      // Case 1: insertions at the same position.
      new_var.add_event(e);
    }
    else if (e.pos == new_var.pos && e.is_deletion() && new_var.is_deletion())
    {
      // Case 2: Deletions at the same position.
      assert(e.ref.size() > new_var.seqs[0].size());

      /// Extend all sequences.
      auto begin_extend = e.ref.cbegin() + new_var.seqs[0].size();

      for (auto & s : new_var.seqs)
        s.append(begin_extend, e.ref.cend());
      ///

      new_var.add_event(e);
    }
    else if (e.pos == new_var.pos && e.is_snp() && new_var.is_snp())
    {
      // Case 3: both are SNPs at the same position.
      new_var.add_event(e);
    }
    else
    {
      // Otherwise, they are two different events.
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


void
Skyr::populate_variants_with_calls()
{
  for (long i = 0; i < static_cast<long>(edits.size()); ++i)
  {
    assert(edits.size() == seqs.size());
    auto const & seq_edits = edits[i];
    uint32_t del_reach = 0;

    for (auto & var : vars)
    {
      assert(var.seqs.size() >= 2);

      auto get_call =
        [del_reach](Variant & variant,
                    std::set<Event2> const & seq_edits)
        {
          std::size_t call = 0;               // Call reference by default
          auto find_it = std::find_if(variant.event_to_allele.begin(),
                                      variant.event_to_allele.end(),
                                      [&](std::pair<Event2, uint32_t> const & e) -> bool
            {
              return seq_edits.count(e.first);
            });

          // Check if the sequence had any event at this variant location
          if (find_it != variant.event_to_allele.end())
          {
            call = find_it->second;
          }
          else if (variant.pos < del_reach)
          {
            assert(variant.seqs.back() == "*");
            // Call asterisk if the variant is deleted by a previous deletion
            //if (variant.seqs[variant.seqs.size() - 1] != "*")
            //  variant.seqs.push_back("*"); // Add asterisk allele if we need to

            call = variant.seqs.size() - 1;
          }

          return call;
        };

      uint16_t const call = get_call(var, seq_edits);
      var.add_call(call);

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
  }
}


} // namespace SIMDPP_ARCH_NAMESPACE
} // namespace paw

#endif // IMPLEMENT_PAW
