#pragma once
#include "rc_man.h"


namespace tinker {
/**
 * \page kmath  Mathematical Algorithm Keywords
 *
 * ### RANDOMSEED [integer]
 * Followed by an integer value, this keyword sets the initial seed value for
 * the random number generator used by Tinker. Setting RANDOMSEED to the same
 * value as an earlier run will allow exact reproduction of the earlier
 * calculation. (Note that this will not hold across different machine types.)
 * RANDOMSEED should be set to a positive integer less than about 2 billion. In
 * the absence of the RANDOMSEED keyword the seed is chosen "randomly" based
 * upon the number of seconds that have elapsed in the current decade.
 */


/**
 * \ingroup rand
 * \brief Returns a random number on `[0,1]` from a uniform distribution.
 */
template <class T>
T random();


/**
 * \ingroup rand
 * \brief Returns a random number from a normal Gaussian distribution
 * with a mean of zero and a variance of one.
 */
template <class T>
T normal();


/**
 * \ingroup rand
 * \brief Returns a random number as if the square of `k` independent standard
 * normal `N(0,1)` random variables were aggregated.
 *
 * \note There is a bug in some (not so) early
 * `g++` `random/chi_squared_distribution` implementations
 * [<a href="https://stackoverflow.com/questions/48248565">ref 1</a> and
 * <a href="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=83833">ref 2</a>],
 * thus gamma distribution is used here via the following equations:
 *
 * \f[ p_{\chi^2}(x|n) = \frac{x^{n/2-1} exp(-x/2)}{2^{n/2} \Gamma(n/2)} \f]
 * \f[ p_\Gamma(x|a,b) = \frac{x^{a-1}   exp(-x/b)}{b^a     \Gamma(a)}   \f]
 * \f[ p_{\chi^2}(n) = p_\Gamma(n/2,2)                                   \f]
 */
template <class T = real>
T chi_squared(int k);


void random_data(rc_op);
}
