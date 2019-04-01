#include <vector>
#include <cmath>

#include "logdouble.h"

logdouble sum(const std::vector<logdouble> &e) {
  // nothing to do if empty
  if (e.empty())
    return 0;

  // find the maximum exponent
  std::vector<logdouble>::const_iterator i = e.begin();
  double l = i->l;
  while (++i != e.end()) {
    if (i->l > l)
      l = i->l;
  }
  if (std::isinf(l) && l < 0)
    return 0; // one of the terms is zero

  double p = 0;
  for (std::vector<logdouble>::const_iterator i = e.begin(); i != e.end(); ++i) {
    p += exp(i->l - l) * i->p;
  }
  return logdouble(l, p);
}
