#pragma once

#include <array>
#include <memory>

#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>
#include <paw/align/alignment_cache.hpp>

#include <simdpp/simd.h>


namespace paw
{

// NOTE: AlignmentResults contain stuff which is irrelevant of type of SIMD instructions used, therefore
// it is (and should not be) in the SIMDPP_ARCH_NAMESPACE namespace
template <typename Tuint>
struct AlignmentResults
{
  long score{0};
  long query_begin{0};
  long query_end{0};
  long database_begin{0};
  long database_end{0};
  long clip_begin{0};
  long clip_end{0};
  std::unique_ptr<std::pair<std::string, std::string> > aligned_strings_ptr;

public:
  template <typename Tseq>
  inline void
  get_aligned_strings(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
                      Tseq const & q,
                      Tseq const & d);

  std::pair<long, long> inline
  get_database_begin_end(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache) const;

  template <typename Tseq>
  std::pair<long, long> inline
  apply_clipping(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
                 Tseq const & q,
                 Tseq const & d,
                 Tuint match,
                 Tuint mismatch,
                 Tuint gap_open,
                 Tuint gap_extend,
                 Tuint clip);


  inline void clear();
  inline void reset(paw::SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> * cache);
};


template <typename Tuint>
template <typename Tseq>
inline void
AlignmentResults<Tuint>::get_aligned_strings(paw::SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
                                             Tseq const & q,
                                             Tseq const & d)
{
  long i = database_end;
  long j = query_end;

  assert(j <= static_cast<long>(q.size()));
  assert(i <= static_cast<long>(d.size()));

  std::pair<std::string, std::string> s;
  //std::cerr << database_end << " " << query_end << " " << d.size() << " " << q.size() << "\n";
  /////*
  //if (query_end < static_cast<long>(q.size()))
  //{
  //  s.first.append(std::string(q.rbegin(), q.rbegin() + q.size() - query_end));
  //  s.second.append(d.size() - query_end, '-');
  //
  //  //del_count = 1;
  //  //del_ext_count = d.size() - query_end - 1;
  //}
  ////*/
  //return s;

  auto add_del = [&]()
                 {
                   assert(j > 0);
                   assert(j <= static_cast<long>(q.size()));
                   //std::cout << "DEL SELECTED " << q[j - 1] << " " << j << "\n";
                   s.first.push_back(q[j - 1]);
                   s.second.push_back('-');
                   --j;
                 };

  auto add_del_ext = [&]()
                     {
                       assert(j > 0);
                       assert(j <= static_cast<long>(q.size()));
                       //std::cout << "DELE SELECTED " << q[j - 1] << " " << j << "\n";
                       s.first.push_back(q[j - 1]);
                       s.second.push_back('-');
                       --j;
                     };

  auto add_ins = [&]()
                 {
                   assert(i > 0);
                   assert(i <= static_cast<long>(d.size()));
                   //std::cout << "INS SELECTED " << d[i - 1] << " " << i << "\n";
                   s.first.push_back('-');
                   s.second.push_back(d[i - 1]);
                   --i;
                 };

  auto add_ins_ext = [&]()
                     {
                       assert(i > 0);
                       assert(i <= static_cast<long>(d.size()));
                       //std::cout << "INSE SELECTED " << d[i - 1] << " " << i << "\n";
                       s.first.push_back('-');
                       s.second.push_back(d[i - 1]);
                       --i;
                     };

  auto add_sub = [&]()
                 {
                   assert(j > 0l);
                   assert(j <= static_cast<long>(q.size()));
                   assert(i > 0l);
                   assert(i <= static_cast<long>(d.size()));
                   //std::cout << "SUB SELECTED " << q[j - 1] << " " << d[i - 1] << "\n";

                   s.first.push_back(q[j - 1]);
                   s.second.push_back(d[i - 1]);
                   --i;
                   --j;
                 };

  while (i > 0 || j > 0)
  {
    assert(s.first.size() == s.second.size());

    if (j == 0)
    {
      while (i > 1)
        add_ins_ext();

      add_ins();
      break;
    }

    assert(i >= 0);
    assert(j >= 0);
    long const v = j % aln_cache.mB.t;
    long const e = j / aln_cache.mB.t;

    if (i == 0)
    {
      assert(j > 0);
      add_del();

      while (j > 0)
        add_del_ext();
    }
    else if (aln_cache.mB.is_del(i - 1, v, e))
    {
      while (j > 1 && aln_cache.mB.is_del_extend(i - 1, j % aln_cache.mB.t, j / aln_cache.mB.t))
        add_del_ext();

      assert(j > 0);
      add_del();
    }
    else if (aln_cache.mB.is_ins(i - 1, v, e))
    {
      while (i > 1 && aln_cache.mB.is_ins_extend(i - 1, v, e))
        add_ins_ext();

      assert(i > 0);
      add_ins();
    }
    else
    {
      add_sub();
    }
  }

  aligned_strings_ptr = std::unique_ptr<std::pair<std::string, std::string> >(
    new std::pair<std::string, std::string>(
      {std::string(s.first.rbegin(), s.first.rend()),
       std::string(s.second.rbegin(), s.second.rend())
      }
      ));
}


template <typename Tuint>
std::pair<long, long> inline
AlignmentResults<Tuint>::get_database_begin_end(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache) const
{
  long i = database_end;
  long j = query_end;

  std::pair<long, long> res = {0, database_end};

  while (i > 0 || j > 0)
  {
    if (j == 0)
    {
      res.first = i;
      break;
    }

    assert(i >= 0);
    assert(j >= 0);
    long const v = j % aln_cache.mB.t;
    long const e = j / aln_cache.mB.t;

    if (i == 0)
    {
      assert(j > 0);
      --j;

      while (j > 0)
        --j;
    }
    else if (aln_cache.mB.is_del(i - 1, v, e))
    {
      while (j > 1 && aln_cache.mB.is_del_extend(i - 1, j % aln_cache.mB.t, j / aln_cache.mB.t))
        --j;

      assert(j > 0);
      --j;
    }
    else if (aln_cache.mB.is_ins(i - 1, v, e))
    {
      while (i > 1 && aln_cache.mB.is_ins_extend(i - 1, v, e))
        --i;

      assert(i > 0);
      --i;

      if (j == query_end)
        res.second = i;
    }
    else
    {
      --i;
      --j;
    }
  }

  return res;
}


template <typename Tuint>
template <typename Tseq>
std::pair<long, long> inline
AlignmentResults<Tuint>::apply_clipping(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
                                        Tseq const & q,
                                        Tseq const & d,
                                        Tuint match,
                                        Tuint mismatch,
                                        Tuint gap_open,
                                        Tuint gap_extend,
                                        Tuint clip)
{
  long i = database_end;
  long j = query_end;

#ifndef NDEBUG
  long const old_score = score;
#endif // NDEBUG

  long tmp_score{0l};
  long best_begin_clip_improvement{0};

  assert(j == static_cast<long>(q.size()));
  assert(i == static_cast<long>(d.size()));

  std::pair<long, long> res = {0, query_end};

  while (i > 0 || j > 0)
  {
    if (j == 0)
    {
      //res.first = i;
      break;
    }

    assert(i >= 0);
    assert(j >= 0);
    long const v = j % aln_cache.mB.t;
    long const e = j / aln_cache.mB.t;

    if (i == 0)
    {
      assert(j > 0);
      --j;

      while (j > 0)
        --j;
    }
    else if (aln_cache.mB.is_del(i - 1, v, e))
    {
      while (j > 1 && aln_cache.mB.is_del_extend(i - 1, j % aln_cache.mB.t, j / aln_cache.mB.t))
      {
        tmp_score -= static_cast<long>(gap_extend);
        --j;
      }

      assert(j > 0);
      tmp_score -= static_cast<long>(gap_open);
      --j;
    }
    else if (aln_cache.mB.is_ins(i - 1, v, e))
    {
      while (i > 1 && aln_cache.mB.is_ins_extend(i - 1, v, e))
      {
        if (j < query_end)
          tmp_score -= static_cast<long>(gap_extend);

        --i;
      }

      assert(i > 0);
      --i;

      if (j < query_end)
        tmp_score -= static_cast<long>(gap_open);
    }
    else
    {
      --i;
      --j;

      if (q[j] == d[i] || q[j] == 'N' || d[i] == 'N')
      {
        // Check clip of end
        if (tmp_score < 0 - static_cast<long>(clip))
        {
          res.second = j + 1;
          score += 0 - static_cast<long>(clip) - tmp_score;
          //std::cout << "BETTER END CLIP " << tmp_score << " > " << old_score << " " <<  i << "," << j << "\n";
          tmp_score = 0l - static_cast<long>(clip);
        }

        // giff match
        tmp_score += static_cast<long>(match);

        // Check clip of begin
        if (tmp_score - static_cast<long>(clip) > score)
        {
          long const diff = tmp_score - static_cast<long>(clip) - score;

          if (diff > best_begin_clip_improvement)
          {
            best_begin_clip_improvement = diff;
            res.first = j;
          }

          //std::cout << "BETTER BEGIN CLIP " << score << " > " << score << " diff=" << diff << "\n";
        }
      }
      else
      {
        //std::cout << "mismatch " << i << "," << j << "\n";
        tmp_score -= static_cast<long>(mismatch);
      }
    }
  }

  //std::cout << "new_score, old_score, tmp_score = " << score << " " << old_score << " " << tmp_score << "\n";
  assert(tmp_score >= old_score);
  score = tmp_score + best_begin_clip_improvement;
  return res;
}


#if defined(IMPLEMENT_PAW)

template <typename Tuint>
void inline
AlignmentResults<Tuint>::clear()
{
  //aln_cache.mB = Backtrack<Tuint>();
  score = 0;
  query_end = 0;
  database_end = 0;
}


#endif // defined(IMPLEMENT_PAW)


} // namespace paw
