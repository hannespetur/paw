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


struct g
{
  g(float seed = 0.0f, float step = 1.0f)
    : i_(seed)
    , step_(step)
  {}

  template < typename T> T operator()(bs::as_<T>)
  {
    auto z = bs::enumerate<T>(i_, step_);
    i_+= bs::cardinal_of<T>::value*step_;
    return z;
  }

  float i_;
  float step_;
};


struct gstd
{
  gstd(float seed = 0.0f, float step = 1.0f)
    : i_(seed)
    ,  step_(step)
  {}

  float operator()()
  {
    float z = i_;
    i_+= step_;
    return z;
  }

  float i_;
  float step_;
};
