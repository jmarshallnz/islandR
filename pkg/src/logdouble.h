#pragma once

class logdouble {
  double p; // current value after log offset
  double l; // logarithmic offset

  inline void check_range() {
    const double too_large = 1.e64;
    const double too_small = 1.e-64;
    if (p > too_large || p < too_small) {
      // take out the logarithm to maintain accuraccy
      l += std::log(p);
      p = 1;
    }
  }

public:
  logdouble() : p(1.0), l(0.0) {}
  logdouble(const double &v) : p(v), l(0.0) { check_range(); }
  logdouble(const double &e, const double &v) : p(v), l(e) { check_range(); }

  logdouble& operator *=(const double &v) {
    p *= v;
    check_range();
    return *this;
  }

  logdouble& operator *=(const logdouble &v) {
    l += v.l; p *= v.p;
    check_range();
    return *this;
  }

  logdouble operator *(const double &v) const {
    return logdouble(l, p * v);
  }

  logdouble operator *(const logdouble &v) const {
    return logdouble(l + v.l, p + v.p);
  }

  logdouble& operator /=(const double &v) {
    p /= v;
    check_range();
    return *this;
  }

  logdouble& operator /=(const logdouble &v) {
    l -= v.l; p /= v.p;
    check_range();
    return *this;
  }

  logdouble operator /(const double &v) const {
    return logdouble(l, p / v);
  }

  logdouble operator /(const logdouble &v) const {
    return logdouble(l - v.l, p / v.p);
  }

  logdouble& operator ^=(const double &v) {
    l = log() * v;
    p = 1;
    return *this;
  }

  logdouble operator ^(const double &v) const {
    return logdouble(log() * v, 1);
  }

  double log() const {
    return l + std::log(p);
  }

  friend logdouble sum(const std::vector<logdouble> &e);
};
