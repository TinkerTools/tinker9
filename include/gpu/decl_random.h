#ifndef TINKER_GPU_DECL_RANDOM_H_
#define TINKER_GPU_DECL_RANDOM_H_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void random_data(rc_t rc);

/// @return  a random number on [0,1] from a uniform distribution
double random_double();
float random_float();
template <class T = real>
T random() {
  if_constexpr(std::is_same<T, double>::value) { random_double(); }
  else if_constexpr(std::is_same<T, float>::value) {
    random_float();
  }
  else {
    assert(false);
  }
}

/**
 * @return
 * a random number from a normal Gaussian distribution with a mean of zero and a
 * variance of one
 */
double normal_double();
float normal_float();
template <class T = real>
T normal() {
  if_constexpr(std::is_same<T, double>::value) { normal_double(); }
  else if_constexpr(std::is_same<T, float>::value) {
    normal_float();
  }
  else {
    assert(false);
  }
}

/**
 * @return
 * a random number as if the square of k independent standard normal (N(0,1))
 * random variables were aggregated; notice that there is a bug in some (not so)
 * early g++ random/chi_squared_distribution implementations, gamma distribution
 * is used instead via the following definitions:
 *
 * chi2_dist  p(x|n)   = x^(n/2-1) exp(-x/2) /GammaFunc(n/2) /2^(n/2)
 *
 * gamma_dist p(x|a,b) = x^(a-1)   exp(-x/b) /GammaFunc(a)   /b^a
 *
 * chi2_dist(n) = gamma_dist(n/2,2)
 */
double chi_squared_double(int k);
float chi_squared_float(int k);
template <class T = real>
T chi_squared(int k) {
  if_constexpr(std::is_same<T, double>::value) { chi_squared_double(k); }
  else if_constexpr(std::is_same<T, float>::value) {
    chi_squared_float(k);
  }
  else {
    assert(false);
  }
}
}
TINKER_NAMESPACE_END

#endif
