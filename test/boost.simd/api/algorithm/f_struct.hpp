#pragma once

template<typename T>
struct f_ge
{
  using p_t = pack<T>;

  f_ge(const T & n)
      : n_(n)
      , vn_(n)
  {}

  bool                  operator()(T a)           const {return a >= n_;  }
  bs::as_logical_t<p_t> operator()(const p_t & a) const {return a >= vn_; }
  T n_;
  p_t vn_;
};


template<typename T>
struct f_lt
{
  using p_t = pack<T>;

  f_lt(const T & n)
      : n_(n)
      , vn_(n)
  {}

  bool                  operator()(T a)           const {return a < n_;  }
  bs::as_logical_t<p_t> operator()(const p_t & a) const {return a < vn_; }
  T n_;
  p_t vn_;
};
