#pragma once

#include <cstdint> // uint8_t
#include <set> // std::set<T>
#include <string> // std::string
#include <vector> // std::vector<T>

#include <paw/align/pairwise_alignment.hpp>
#include <paw/align/event.hpp>
#include <paw/align/variant.hpp>


namespace paw
{

class Skyr
{
public:
  std::vector<std::string> seqs;
  std::vector<std::set<Event2> > edits;   // Edits found for each alignment
  std::multiset<Event2> all_edits;
  std::vector<Variant> vars;

  explicit Skyr(std::vector<std::string> const & _seqs);
  explicit Skyr(std::vector<std::vector<char> > const & _seqs);

  void find_all_edits(bool const is_normalize);
  void find_variants_from_edits();
  void populate_variants_with_calls(bool use_asterisks = true);
  Variant merge_variants(long start_v, long end_v);
  std::vector<Variant> split_variants(long const T, std::string const & pad_base = "");
};


} // namespace paw


#if defined(IMPLEMENT_PAW) || defined(__JETBRAINS_IDE__)


#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>


namespace
{

long
_get_variant_call(paw::Variant & variant,
                  std::set<paw::Event2> const & seq_edits,
                  long del_reach,
                  bool use_asterisks = true)
{
  long call = 0; // Call reference by default
  auto find_it = std::find_if(variant.event_to_allele.begin(),
                              variant.event_to_allele.end(),
                              [&](std::pair<paw::Event2, uint32_t> const & e) -> bool
    {
      return seq_edits.count(e.first);
    });

  // Check if the sequence had any event at this variant location
  if (find_it != variant.event_to_allele.end())
  {
    call = find_it->second;
  }
  else if (use_asterisks && variant.pos < del_reach)
  {
    // Call asterisk if the variant is deleted by a previous deletion
    if (variant.seqs[variant.seqs.size() - 1] != "*")
      variant.seqs.push_back("*"); // Add asterisk allele if we need to

    call = variant.seqs.size() - 1;
  }

  return call;
}


long
_const_get_variant_call(paw::Variant const & variant,
                        std::set<paw::Event2> const & seq_edits,
                        long del_reach)
{
  long call = 0; // Call reference by default
  auto find_it = std::find_if(variant.event_to_allele.begin(),
                              variant.event_to_allele.end(),
                              [&](std::pair<paw::Event2, uint32_t> const & e) -> bool
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
    call = -1;
  }

  return call;
}


} // anon namespace


namespace paw
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


void
Skyr::find_all_edits(bool const is_normalize)
{
  using Tuint = uint16_t;
  AlignmentOptions<Tuint> opts;
  opts.get_aligned_strings = true;
  opts.set_match(1).set_mismatch(4).set_gap_open(7).set_gap_extend(1);
  all_edits.clear();

  for (int i{1}; i < static_cast<int>(seqs.size()); ++i)
  {
    assert(i < static_cast<int>(edits.size()));

    paw::pairwise_alignment<std::string, Tuint>(seqs[0], seqs[i], opts);
    auto ar = opts.get_alignment_results();
    assert(ar);

    assert(ar->aligned_strings_ptr);
    auto const & aligned_strings = *(ar->aligned_strings_ptr);

    // std::cout << aligned_strings.first << '\n' << aligned_strings.second << '\n';
    edits[i] = get_edit_script(aligned_strings, is_normalize, false);
    // for (auto it = edits[i].begin(); it != edits[i].end(); ++it)
    // {
    //   std::cout << it->pos << " " << it->ref << " " << it->alt << '\n';
    // }
    all_edits.insert(edits[i].begin(), edits[i].end());
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
        continue; // Same as next edit, we can skip it
    }

    auto const & e = *it;

    if (!new_var.has_sequences())
    {
      new_var = {static_cast<uint32_t>(e.pos), {e.ref}};
      new_var.add_event(e);
      continue;
    }

    if (e.pos == new_var.pos)
    {
      if (e.is_indel() && new_var.is_indel())
      {
        // Case 1: indel at the same position.
        assert(e.ref.size() > new_var.seqs[0].size());

        // Extend all variant sequences if needed
        {
          auto const begin_extend = std::next(e.ref.cbegin(), new_var.seqs[0].size());
          auto const end_extend = e.ref.cend();

          for (auto & s : new_var.seqs)
            s.append(begin_extend, end_extend);
        }

        new_var.add_event(e);
      }
      else if (e.is_snp() && new_var.is_snp())
      {
        // Case 2: both are SNPs at the same position.
        new_var.add_event(e);
      }
      else
      {
        // Otherwise, they are two different events.
        vars.push_back(std::move(new_var)); // push old variant
        new_var = {static_cast<uint32_t>(e.pos), {e.ref} };
        new_var.add_event(e);
      }
    }
    else
    {
      // Otherwise, they are two different events.
      vars.push_back(std::move(new_var)); // push old variant
      new_var = {static_cast<uint32_t>(e.pos), {e.ref} };
      new_var.add_event(e);
    }
  }

  if (new_var.seqs.size() > 1)
  {
    vars.push_back(std::move(new_var));
  }
}


void
Skyr::populate_variants_with_calls(bool const use_asterisks)
{
  for (long i = 0; i < static_cast<long>(edits.size()); ++i)
  {
    assert(edits.size() == seqs.size());
    auto const & seq_edits = edits[i];
    long del_reach = -1;

    for (auto & var : vars)
    {
      assert(var.seqs.size() >= 2);
      long const call = _get_variant_call(var, seq_edits, del_reach, use_asterisks);
      var.add_call(call);

      if (var.is_deletion())
      {
        assert(var.seqs[0].size() >= var.seqs[call].size());
        del_reach = std::max(del_reach, static_cast<long>(var.pos +
                                                          var.seqs[0].size() -
                                                          var.seqs[call].size()));
      }
    }
  }
}


Variant
Skyr::merge_variants(long v1, long v2)
{
  assert(seqs.size() > 0);
  assert(v1 <= v2);

  //std::string const & ref = seqs[0];
  std::cerr << "split [" << v1 << "," << v2 << "]\n";

  //if (v1 == v2)
  //  return vars[v1];

  Variant new_var;
  new_var.pos = vars[v1].pos;

  // add reference
  {
    long s1 = vars[v1].pos;
    long s2;

    if (v2 + 1l == static_cast<long>(vars.size()))
      s2 = std::string::npos;
    else
      s2 = vars[v2 + 1].pos - s1;

    new_var.seqs.push_back(seqs[0].substr(s1, s2));
  }

  for (long i = 1; i < static_cast<long>(edits.size()); ++i)
  {
    assert(edits.size() == seqs.size());
    auto const & edit = edits[i];

    long shift = 0;
    long s1 = -1;
    long s2 = -1;

    for (auto const & e : edit)
    {
      if ((v2 + 1l) < static_cast<long>(vars.size()) && e.pos >= vars[v2 + 1].pos)
      {
        assert(e.pos == vars[v2 + 1].pos);
        s2 = vars[v2 + 1].pos - s1 + shift;
        std::cerr << "DONE " << e.pos << " " << e.ref << " " << e.alt << std::endl;
        break;
      }

      if (s1 == -1 && e.pos == vars[v1].pos)
      {
        s1 = vars[v1].pos + shift;
      }

      assert(e.alt != "*");
      shift += static_cast<long>(e.alt.size()) - static_cast<long>(e.ref.size());
      std::cerr << "GO " << e.pos << " " << e.ref << " " << e.alt << "\n";
    }

    std::cerr << "info: " << s1 << "," << s2 << "," << shift << std::endl;
    assert(s1 != -1);

    if (s2 == -1)
      new_var.seqs.push_back(seqs[i].substr(s1));
    else
      new_var.seqs.push_back(seqs[i].substr(s1, s2));

    std::cerr << new_var.seqs.back() << " (" << i << ")" << std::endl;
  }

  return new_var;
}


std::vector<Variant>
Skyr::split_variants(long const T, std::string const & pad_base)
{
  std::vector<Variant> new_variants;

  if (vars.size() == 0)
    return new_variants; // No variants, nothing to do

  // Process reference
  {
    long var_reach = -1;//vars[0].get_max_reach();
    long ref_s = 0; // string index
    std::string original_seq = seqs[0];

    for (long v = 0; v < static_cast<long>(vars.size()); ++v)
    {
      auto const & var = vars[v];
      assert(var.seqs.size() >= 2);
      var_reach = std::max(var_reach, var.get_max_reach());

      // Check if there are more than T base till the next variant
      if (v + 1l < static_cast<long>(vars.size()))
      {
        if (((var_reach + T) <= vars[v + 1].pos) ||
            ((var_reach - ref_s) > 60 && (var_reach + (T / 2)) < vars[v + 1].pos) ||
            ((var_reach - ref_s) > 80 && (var_reach + (T / 8)) < vars[v + 1].pos) ||
            ((var_reach - ref_s) > 100 && var_reach < vars[v + 1].pos))
        {
          Variant new_var;
          new_var.pos = ref_s;
          new_var.seqs.push_back(original_seq.substr(ref_s, var_reach - ref_s)); //
          new_var.add_call(0);
          new_variants.push_back(std::move(new_var));
          ref_s = vars[v + 1].pos;
        }
      }
    }

    Variant new_var;
    new_var.pos = ref_s;
    new_var.seqs.push_back(original_seq.substr(ref_s));
    new_var.add_call(0);
    new_variants.push_back(std::move(new_var));
  }

  // Process alts
  for (long i = 1; i < static_cast<long>(edits.size()); ++i)
  {
    auto const & edit = edits[i];
    long del_reach = -1;//vars[0].get_max_del_reach();
    long var_reach = -1;//vars[0].get_max_reach();
    long shift = 0;
    long ref_s = 0;
    long s = 0; // string index
    long new_var_index = 0;
    std::string original_seq = seqs[0]; // Reconstructed sequence
    std::string gapped_seq = seqs[0]; // Gapped reconstructed sequence

    for (long v = 0; v < static_cast<long>(vars.size()); ++v)
    {
      auto const & var = vars[v];
      assert(var.seqs.size() >= 2);
      //long const call = _get_variant_call(var, edit, del_reach, true);
      long const call = _const_get_variant_call(var, edit, del_reach);
      var_reach = std::max(var_reach, var.get_max_reach());

      // TODO (maybe)
      // I need to change var.get_max_del_reach() to static_cast<long>(var.pos + var.seqs[0].size() - var.seqs[call].size())
      if (var.is_deletion())
        del_reach = std::max(del_reach, var.get_max_del_reach());

      if (call > 0)
      {
        original_seq = original_seq.replace(var.pos + shift, var.seqs[0].size(), var.seqs[call]);
        shift += var.seqs[call].size() - var.seqs[0].size();
      }

      // Check if there are more than T base till the next variant
      if (v + 1l < static_cast<long>(vars.size()))
      {
        if (((var_reach + T) <= vars[v + 1].pos) ||
            ((var_reach - ref_s) > 60 && (var_reach + (T / 2)) < vars[v + 1].pos) ||
            ((var_reach - ref_s) > 80 && (var_reach + (T / 8)) < vars[v + 1].pos) ||
            ((var_reach - ref_s) > 100 && var_reach < vars[v + 1].pos))
        {
          assert(new_var_index < static_cast<long>(new_variants.size()));
          std::string new_seq = original_seq.substr(s, var_reach + shift - s);
          paw::Variant & new_var = new_variants[new_var_index];

          {
            auto & new_seqs = new_var.seqs;

            // Only add sequence if we did not see it before
            auto find_it = std::find(new_seqs.begin(), new_seqs.end(), new_seq);
            new_var.add_call(std::distance(new_seqs.begin(), find_it));

            if (find_it == new_seqs.end())
            {
              new_var.seqs.push_back(std::move(new_seq)); // new sequence
            }
          }

          ++new_var_index;
          s = vars[v + 1].pos + shift;
          ref_s = vars[v + 1].pos;
        }
      }

      //if (v + 1l < static_cast<long>(vars.size()) && var_reach + T <= vars[v + 1].pos)
      //{
      //
      //}
    }

    std::string new_seq = original_seq.substr(s);
    assert(new_var_index < static_cast<long>(new_variants.size()));
    paw::Variant & new_var = new_variants[new_var_index];

    {
      auto & new_seqs = new_var.seqs;

      // Only add sequence if we did not see it before
      auto find_it = std::find(new_seqs.begin(), new_seqs.end(), new_seq);
      new_var.add_call(std::distance(new_seqs.begin(), find_it));

      if (find_it == new_seqs.end())
      {
        new_var.seqs.push_back(std::move(new_seq)); // new sequence
      }
    }
  }

  for (auto & new_var : new_variants)
  {
    new_var.add_base_to_front(seqs[0], pad_base);
    new_var.trim_sequences();
  }

  return new_variants;
}


/*
long prev_v = 0;

for (long v = 1; v < static_cast<long>(vars.size()); ++v)
{
  auto & var = vars[v];
  var.print_seqs(std::cerr);

  if (var_reach + T <= var.pos)
  {
    Variant new_var = merge_variants(prev_v, v - 1);
    prev_v = v;
    new_variants.push_back(std::move(new_var));
    }*/

//assert(var.seqs.size() >= 2);
//long const call = _get_variant_call(var, edit, del_reach, true); // use_asterisks is true
//var.add_call(call);
//
//var.print_seqs(std::cerr);
//std::cerr << "call = " << var.seqs[call] << "\n";

/* // Check if we need to replace strings in the reconstructed sequence
if (call != 0l && var.seqs[call] != "*")
{
  long new_pos = var.pos + shift;
  original_seq = original_seq.replace(new_pos, var.seqs[0].size(), var.seqs[call]);
  shift += var.seqs[call].size() - var.seqs[0].size();
}

if (var.seqs[call] != "*")
{
  long new_gapped_pos = var.pos + gapped_shift;
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
*/

/*
    var_reach = std::max(var_reach, var.get_max_reach());

    if (var.is_deletion())
      del_reach = std::max(del_reach, var.get_max_del_reach());

    if (prev_v < static_cast<long>(vars.size()))
      vars[prev_v].print_seqs(std::cerr);
  }

  {
    long v = static_cast<long>(vars.size());

    assert(prev_v < v);
    Variant new_var = merge_variants(prev_v, v - 1);
    new_variants.push_back(std::move(new_var));
  }

  for (auto const & new_var: new_variants)
  {
    new_var.print_seqs(std::cerr);
  }
*/

} // namespace paw


#endif // IMPLEMENT_PAW
