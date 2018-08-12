#pragma once

#include <paw/align/libsimdpp_backtracker.hpp>
#include <paw/align/libsimdpp_utils.hpp>

#include <simdpp/simd.h>


namespace paw
{

template <typename Tuint>
struct AlignmentResults
{
  long score{0};
  SIMDPP_ARCH_NAMESPACE::Backtrack<Tuint> mB;
  long query_end{0};
  long database_end{0};
  long t{0};

  AlignmentResults() = default;

  template<typename Tseq>
  std::pair<std::string, std::string> get_aligned_strings(Tseq const & q, Tseq const & d) const;
};


//#if defined(IMPLEMENT_PAW)

template<typename Tuint>
template<typename Tseq>
std::pair<std::string, std::string>
AlignmentResults<Tuint>::get_aligned_strings(Tseq const & q, Tseq const & d) const
{
  long i = database_end;
  long j = query_end;

  assert(j <= (long)q.size());
  assert(i <= (long)d.size());

  std::pair<std::string, std::string> s;

  /*
  if (m < static_cast<long>(q.size()))
  {
    s.first.append(std::string(q.rbegin(), q.rbegin() + q.size() - query_end));
    s.second.append(d.size() - query_end, '-');

    del_count = 1;
    del_ext_count = d.size() - query_end - 1;
  }
  */

  assert(s.first.size() == s.second.size());

  auto add_del = [&]()
  {
    assert(j > 0);
    assert(j <= (long)q.size());
    //std::cout << "DEL SELECTED " << q[j - 1] << " " << j << "\n";
    s.first.push_back(q[j - 1]);
    s.second.push_back('-');
    --j;
  };

  auto add_del_ext = [&]()
  {
    assert(j > 0);
    assert(j <= (long)q.size());
    //std::cout << "DELE SELECTED " << q[j - 1] << " " << j << "\n";
    s.first.push_back(q[j - 1]);
    s.second.push_back('-');
    --j;
  };

  auto add_ins = [&]()
  {
    assert(i > 0);
    assert(i <= (long)d.size());
    //std::cout << "INS SELECTED " << d[i - 1] << " " << i << "\n";
    s.first.push_back('-');
    s.second.push_back(d[i - 1]);
    --i;
  };

  auto add_ins_ext = [&]()
  {
    assert(i > 0);
    assert(i <= (long)d.size());
    //std::cout << "INSE SELECTED " << d[i - 1] << " " << i << "\n";
    s.first.push_back('-');
    s.second.push_back(d[i - 1]);
    --i;
  };

  auto add_sub = [&]()
  {
    assert(j > 0l);
    assert(j <= (long)q.size());
    assert(i > 0l);
    assert(i <= (long)d.size());
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
    long const v = j % t;
    long const e = j / t;

    if (i == 0)
    {
      assert(j > 0);
      add_del();

      while (j > 0)
        add_del_ext();
    }
    else if (mB.is_del(i - 1, v, e))
    {
      while (j > 1 && mB.is_del_extend(i - 1, j % t, j / t))
        add_del_ext();

      assert(j > 0);
      add_del();
    }
    else if (mB.is_ins(i - 1, v, e))
    {
      while (i > 1 && mB.is_ins_extend(i - 1, v, e))
        add_ins_ext();

      assert(i > 0);
      add_ins();
    }
    else
    {
      add_sub();
    }
  }

  std::pair<std::string, std::string> out =
  {
    std::string(s.first.rbegin(), s.first.rend()),
    std::string(s.second.rbegin(), s.second.rend())
  };

  assert(out.first.size() == out.second.size());
  return out;
}

//#endif // defined(IMPLEMENT_PAW)


} // namespace paw