#pragma once

#include <vector>

#include <simdpp/simd.h>


namespace paw
{
namespace SIMDPP_ARCH_NAMESPACE
{


template<typename Tpack>
void
print_pack(Tpack const & pack)
{
  using T = typename Tpack::uint_element_type;

  // Guard for when the pack is empty
  if (pack.length == 0)
  {
    std::cout << "()";
    return;
  }

  std::vector<T, simdpp::aligned_allocator<T, sizeof(T)> > vector;
  vector.resize(pack.length);
  simdpp::store(&vector[0], pack);

  std::cout << "(" << static_cast<long>(vector[0]);

  for (long i = 1; i < static_cast<long>(vector.size()); ++i)
  {
    std::cout << "," << static_cast<long>(vector[i]);
  }

  std::cout << ")" << std::endl;
}


template <typename Trow>
void
print_score_vector_standard(Trow const & vX)
{
  using Tuint = typename Trow::Tuint;
  using Vec = std::vector<Tuint, simdpp::aligned_allocator<Tuint, sizeof(Tuint)> >;
  using Mat = std::vector<Vec>;

  long const t = vX.vectors.size();

  if (t == 0)
    return;

  Vec vec(vX.vectors[0].length, 0);
  Mat m(vX.vectors.size(), vec);

  for (long v = 0; v < t; ++v)
  {
    simdpp::store(&m[v][0], vX.vectors[v]);
  }

  for (std::size_t j = 0; j < vX.n_elements; ++j)
  {
    std::size_t const v = j % t;
    std::size_t const e = j / t;
    std::cout << std::setw(4) << static_cast<uint64_t>(m[v][e]) << " ";
  }

  std::cout << "\n";
}


template <typename Tpack>
void
print_score_vectors(Tpack const & vH,
                    Tpack const & vH_up,
                    Tpack const & vE,
                    Tpack const & vF,
                    Tpack const & vF_up,
                    Tpack const & vW
                    )
{
  std::cout << "Standard H_up  : "; print_score_vector_standard(vH_up);
  std::cout << "Standard H     : "; print_score_vector_standard(vH);
  std::cout << "Standard E     : "; print_score_vector_standard(vE);
  std::cout << "Standard F_up  : "; print_score_vector_standard(vF_up);
  std::cout << "Standard F     : "; print_score_vector_standard(vF);
  std::cout << "Standard W     : "; print_score_vector_standard(vW);
  std::cout << "=====\n";
}


} //namespace SIMDPP_ARCH_NAMESPACE
} // anon namespace
