#ifndef TINKER_RANDOM_H_
#define TINKER_RANDOM_H_

#include "macro.h"
#include <cassert>
#include <type_traits>

TINKER_NAMESPACE_BEGIN
/**
 * @return
 * a random number on [0,1] from a uniform distribution
 */
/// @{
double random_double ();
float random_float ();
template <class T = real>
T random ()
{
   if_constexpr (std::is_same<T, double>::value)
   {
      random_double ();
   }
   else if_constexpr (std::is_same<T, float>::value)
   {
      random_float ();
   }
   else
   {
      assert (false);
   }
}
/// @}

/**
 * @return
 * a random number from a normal Gaussian distribution with a mean of zero and a
 * variance of one
 */
/// @{
double normal_double ();
float normal_float ();
template <class T = real>
T normal ()
{
   if_constexpr (std::is_same<T, double>::value)
   {
      normal_double ();
   }
   else if_constexpr (std::is_same<T, float>::value)
   {
      normal_float ();
   }
   else
   {
      assert (false);
   }
}
/// @}

/**
 * @return
 * a random number as if the square of k independent standard normal N(0,1)
 * random variables were aggregated; notice that there is a bug in some (not so)
 * early @c g++ @c random/chi_squared_distribution implementations
 * [<a href="https://stackoverflow.com/questions/48248565">ref 1</a> and
 * <a href="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=83833">ref 2</a>],
 * gamma distribution is used here via the following equations:
 *
 * @f[ p_{\chi^2}(x|n) = \frac{x^{n/2-1} exp(-x/2)}{2^{n/2} \Gamma(n/2)} @f]
 * @f[ p_\Gamma(x|a,b) = \frac{x^{a-1}   exp(-x/b)}{b^a     \Gamma(a)}   @f]
 * @f[ p_{\chi^2}(n) = p_\Gamma(n/2,2)                                   @f]
 */
/// @{
double chi_squared_double (int k);
float chi_squared_float (int k);
template <class T = real>
T chi_squared (int k)
{
   if_constexpr (std::is_same<T, double>::value)
   {
      chi_squared_double (k);
   }
   else if_constexpr (std::is_same<T, float>::value)
   {
      chi_squared_float (k);
   }
   else
   {
      assert (false);
   }
}
/// @}
TINKER_NAMESPACE_END

#endif
