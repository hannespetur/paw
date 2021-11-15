#pragma once

#include <array>
#include <memory>
#include <limits>

#include <paw/align/alignment_cache.hpp>
#include <paw/align/cigar.hpp>
#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>

#include <simdpp/simd.h>

namespace paw
{
// NOTE: AlignmentResults contain stuff which is irrelevant of type of SIMD instructions used, therefore
// it is (and should not be) in the SIMDPP_ARCH_NAMESPACE namespace
struct AlignmentResults
{
  int32_t score{std::numeric_limits<int32_t>::min()};
  int32_t query_begin{0};
  int32_t query_end{0};
  int32_t database_begin{0};
  int32_t database_end{0};
  int32_t clip_begin{0};
  int32_t clip_end{0};
  std::unique_ptr<std::pair<std::string, std::string>> aligned_strings_ptr{nullptr};
  std::unique_ptr<std::vector<paw::Cigar>> cigar_string_ptr{nullptr};

  AlignmentResults()
  {}

  AlignmentResults(AlignmentResults const & o) = delete;
  AlignmentResults(AlignmentResults && o) noexcept
    : score(o.score)
    , query_begin(o.query_begin)
    , query_end(o.query_end)
    , database_begin(o.database_begin)
    , database_end(o.database_end)
    , clip_begin(o.clip_begin)
    , clip_end(o.clip_end)
    , aligned_strings_ptr(std::forward<std::unique_ptr<std::pair<std::string, std::string>>>(o.aligned_strings_ptr))
    , cigar_string_ptr(std::forward<std::unique_ptr<std::vector<paw::Cigar>>>(o.cigar_string_ptr))
  {
  }

  AlignmentResults& operator=(AlignmentResults const & o) = delete;
  AlignmentResults& operator=(AlignmentResults && o) noexcept
  {
    score = o.score;
    query_begin = o.query_begin;
    query_end = o.query_end;
    database_begin = o.database_begin;
    database_end = o.database_end;
    clip_begin = o.clip_begin;
    clip_end = o.clip_end;
    aligned_strings_ptr = std::move(o.aligned_strings_ptr);
    cigar_string_ptr = std::move(o.cigar_string_ptr);
    return *this;
  }

  ~AlignmentResults() noexcept = default;

public:
  template <typename Tuint>
  inline void get_cigar_string(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache);

  template <typename Tuint, typename Tseq>
  inline void get_aligned_strings(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
                                  Tseq const & q,
                                  Tseq const & d);

  template <typename Tuint>
  std::pair<long, long> inline get_database_begin_end(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache) const;

  template <typename Tuint, typename Tseq>
  std::pair<long, long> inline apply_clipping(SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
                                              Tseq const & q,
                                              Tseq const & d,
                                              Tuint match,
                                              Tuint mismatch,
                                              Tuint gap_open,
                                              Tuint gap_extend,
                                              Tuint clip);

  inline void clear();

  template<typename Tuint>
  inline void reset(paw::SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> * cache);
};

template <typename Tuint>
inline void AlignmentResults::get_cigar_string(paw::SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache)
{
  long i = database_end;
  long j = query_end;

  cigar_string_ptr = std::make_unique<std::vector<paw::Cigar>>();
  std::vector<Cigar> & cigar_string = *cigar_string_ptr;

  if (clip_end < query_end)
  {
    cigar_string.emplace_back(query_end - clip_end, CigarOperation::SOFT_CLIP);
    j = clip_end;
  }

  if (j <= clip_begin)
    return;

  Cigar new_cigar;

  auto add_del = [this](long & j, Cigar & new_cigar, std::vector<Cigar> & cigar_string)
  {
    assert(j > 0);

    if (j > this->clip_end)
    {
      --j;
      return;
    }

    if (new_cigar.operation == CigarOperation::DELETION)
    {
      //std::cout << "Del now: " << new_cigar.count << "\n";
      ++new_cigar.count;
    }
    else
    {
      if (new_cigar.operation != CigarOperation::UNSET)
        cigar_string.push_back(new_cigar);

      new_cigar.operation = CigarOperation::DELETION;
      new_cigar.count = 1;
    }

    --j;
  };

  auto add_del_ext = [this](long & j, Cigar & new_cigar, std::vector<Cigar> & cigar_string)
  {
    assert(j > 0);

    if (j > this->clip_end)
    {
      --j;
      return;
    }

    if (new_cigar.operation == CigarOperation::DELETION)
    {
      ++new_cigar.count;
    }
    else
    {
      if (new_cigar.operation != CigarOperation::UNSET)
        cigar_string.push_back(new_cigar);

      new_cigar.operation = CigarOperation::DELETION;
      new_cigar.count = 1;
    }

    --j;
  };

  auto add_ins = [this](long & i, Cigar & new_cigar, std::vector<Cigar> & cigar_string)
  {
    assert(i > 0);

    if (new_cigar.operation == CigarOperation::INSERTION)
    {
      //std::cout << "Ins now: " << new_cigar.count << "\n";
      ++new_cigar.count;
    }
    else
    {
      if (new_cigar.operation != CigarOperation::UNSET)
        cigar_string.push_back(Cigar(new_cigar));

      new_cigar.operation = CigarOperation::INSERTION;
      new_cigar.count = 1;
    }

    --i;
  };

  auto add_ins_ext = [](long & i, Cigar & new_cigar, std::vector<Cigar> & cigar_string)
  {
    assert(i > 0);

    if (new_cigar.operation == CigarOperation::INSERTION)
    {
      //std::cout << "Ins now: " << new_cigar.count << "\n";
      ++new_cigar.count;
    }
    else
    {
      if (new_cigar.operation != CigarOperation::UNSET)
        cigar_string.push_back(Cigar(new_cigar));

      new_cigar.operation = CigarOperation::INSERTION;
      new_cigar.count = 1;
    }

    --i;
  };

  auto add_sub = [this](long & i, long & j, Cigar & new_cigar, std::vector<Cigar> & cigar_string)
  {
    assert(j > 0l);
    assert(i > 0l);

    if (j > this->clip_end)
    {
      --j;
      return;
    }

    if (new_cigar.operation == CigarOperation::MATCH)
    {
      ++new_cigar.count;
    }
    else
    {
      if (new_cigar.operation != CigarOperation::UNSET)
        cigar_string.push_back(Cigar(new_cigar));

      new_cigar.operation = CigarOperation::MATCH;
      new_cigar.count = 1;
    }

    --i;
    --j;
  };

  while ((i > database_begin || j > clip_begin) && (clip_begin == 0 || j > clip_begin))
  {
    if (j == 0)
    {
      while (i > (database_begin + 1))
        add_ins_ext(i, new_cigar, cigar_string);

      add_ins(i, new_cigar, cigar_string);
      break;
    }

    assert(i >= 0);
    assert(j >= 0);
    long const v = j % aln_cache.mB.t;
    long const e = j / aln_cache.mB.t;

    if (i == database_begin)
    {
      assert(j > clip_begin);
      add_del(j, new_cigar, cigar_string);

      while (j > clip_begin)
        add_del_ext(j, new_cigar, cigar_string);
    }
    else if (aln_cache.mB.is_del(i - 1, v, e))
    {
      while (j > clip_begin + 1 && aln_cache.mB.is_del_extend(i - 1, j % aln_cache.mB.t, j / aln_cache.mB.t))
        add_del_ext(j, new_cigar, cigar_string);

      assert(j > clip_begin);
      add_del(j, new_cigar, cigar_string);
    }
    else if (aln_cache.mB.is_ins(i - 1, v, e))
    {
      if (j > clip_end)
      {
        while (i > (database_begin + 1) && aln_cache.mB.is_ins_extend(i - 1, v, e))
          --i;

        assert(i > database_begin);
        --i;
      }
      else
      {
        while (i > (database_begin + 1) && aln_cache.mB.is_ins_extend(i - 1, v, e))
        {
          if (j > clip_end)
            --i;
          else
            add_ins_ext(i, new_cigar, cigar_string);
        }

        assert(i > database_begin);

        if (j > clip_end)
          --i;
        else
          add_ins(i, new_cigar, cigar_string);
      }

    }
    else
    {
      add_sub(i, j, new_cigar, cigar_string);
    }
  }

  if (new_cigar.operation != CigarOperation::UNSET)
    cigar_string.push_back(new_cigar);

  if (clip_begin > 0)
    cigar_string.emplace_back(clip_begin, CigarOperation::SOFT_CLIP);
}

template <typename Tuint, typename Tseq>
inline void AlignmentResults::get_aligned_strings(paw::SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
                                                         Tseq const & /*q*/,
                                                         Tseq const & d)
{
  long i = database_end;
  long j = query_end;

  auto const & q = aln_cache.query;
  assert(j <= static_cast<long>(q.size()));
  assert(i <= static_cast<long>(d.size()));

  std::pair<std::string, std::string> s;
  // std::cerr << database_end << " " << query_end << " " << d.size() << " " << q.size() << "\n";
  /////*
  // if (query_end < static_cast<long>(q.size()))
  //{
  //  s.first.append(std::string(q.rbegin(), q.rbegin() + q.size() - query_end));
  //  s.second.append(d.size() - query_end, '-');
  //
  //  //del_count = 1;
  //  //del_ext_count = d.size() - query_end - 1;
  //}
  ////*/
  // return s;

  auto add_del = [&]()
  {
    assert(j > 0);
    assert(j <= static_cast<long>(q.size()));
    // std::cout << "DEL SELECTED " << q[j - 1] << " " << j << "\n";
    s.first.push_back(q[j - 1]);
    s.second.push_back('-');
    --j;
  };

  auto add_del_ext = [&]()
  {
    assert(j > 0);
    assert(j <= static_cast<long>(q.size()));
    // std::cout << "DELE SELECTED " << q[j - 1] << " " << j << "\n";
    s.first.push_back(q[j - 1]);
    s.second.push_back('-');
    --j;
  };

  auto add_ins = [&]()
  {
    assert(i > 0);
    assert(i <= static_cast<long>(d.size()));
    // std::cout << "INS SELECTED " << d[i - 1] << " " << i << "\n";
    s.first.push_back('-');
    s.second.push_back(d[i - 1]);
    --i;
  };

  auto add_ins_ext = [&]()
  {
    assert(i > 0);
    assert(i <= static_cast<long>(d.size()));
    // std::cout << "INSE SELECTED " << d[i - 1] << " " << i << "\n";
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
    // std::cout << "SUB SELECTED " << q[j - 1] << " " << d[i - 1] << "\n";

    s.first.push_back(q[j - 1]);
    s.second.push_back(d[i - 1]);
    --i;
    --j;
  };

  while (i > database_begin || j > 0)
  {
    assert(s.first.size() == s.second.size());

    if (j == 0)
    {
      while (i > (database_begin + 1))
        add_ins_ext();

      add_ins();
      break;
    }

    assert(i >= 0);
    assert(j >= 0);
    long const v = j % aln_cache.mB.t;
    long const e = j / aln_cache.mB.t;

    if (i == database_begin)
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
      while (i > (database_begin + 1) && aln_cache.mB.is_ins_extend(i - 1, v, e))
        add_ins_ext();

      assert(i > database_begin);
      add_ins();
    }
    else
    {
      add_sub();
    }
  }

  aligned_strings_ptr = std::unique_ptr<std::pair<std::string, std::string>>(new std::pair<std::string, std::string>(
    {std::string(s.first.rbegin(), s.first.rend()), std::string(s.second.rbegin(), s.second.rend())}));
}

template <typename Tuint>
std::pair<long, long> inline AlignmentResults::get_database_begin_end(
  SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache) const
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

template <typename Tuint, typename Tseq>
inline std::pair<long, long> AlignmentResults::apply_clipping(
  SIMDPP_ARCH_NAMESPACE::AlignmentCache<Tuint> & aln_cache,
  Tseq const & q,
  Tseq const & d,
  Tuint match,
  Tuint mismatch,
  Tuint gap_open,
  Tuint gap_extend,
  Tuint clip)
{
  //#ifndef NDEBUG
  //long const old_score = score;
  long best_begin_clip_improvement{0};
  //#endif // NDEBUG

  long tmp_score{0l};
  std::pair<long, long> res = {0, query_end};

  {
    long i = database_end;
    long j = query_end;

    assert(j == static_cast<long>(q.size()));
    assert(i == static_cast<long>(d.size()));

    while (i > 0 || j > 0)
    {
      if (j == 0)
      {
        // res.first = i;
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
            score -= static_cast<long>(clip) + tmp_score;
            //std::cout << "BETTER END CLIP " << tmp_score << " > " << old_score << " " << i << "," << j << "\n";

            // adjust best_begin_clip based on this new end clip
            best_begin_clip_improvement += static_cast<long>(clip) + tmp_score;

            if (best_begin_clip_improvement <= 0)
            {
              //if (res.first > 0)
              //  std::cout << "Begin clip is now bad: " << best_begin_clip_improvement << " " << res.first << "\n";

              best_begin_clip_improvement = 0;
              res.first = 0;
            }

            tmp_score = -static_cast<long>(clip);
          }

          // giff match
          tmp_score += static_cast<long>(match);

          // Check clip of begin
          if (tmp_score - static_cast<long>(clip) > score)
          {
            long const diff = tmp_score - static_cast<long>(clip) - score;

            if (diff > best_begin_clip_improvement)
            {
              //std::cout << "BETTER BEGIN CLIP " << diff << " > " << best_begin_clip_improvement
              //          << " tmp_score=" << tmp_score << " score=" << score << "\n";

              best_begin_clip_improvement = diff;
              res.first = j;
            }
          }
        }
        else
        {
          tmp_score -= static_cast<long>(mismatch);
        }
      }
    }
  }

  //std::cout << "new_score, old_score, tmp_score = " << score << " " << old_score << " " << tmp_score << "\n";
  //assert(tmp_score >= old_score);
  score = tmp_score + best_begin_clip_improvement;
  return res;
}

#if defined(IMPLEMENT_PAW)

void inline AlignmentResults::clear()
{
  // aln_cache.mB = Backtrack<Tuint>();
  score = 0;
  query_end = 0;
  database_end = 0;
}

#endif // defined(IMPLEMENT_PAW)

} // namespace paw
