#ifndef TINKER_IO_PRINT_H_
#define TINKER_IO_PRINT_H_

#include "macro.h"
#include <ext/fmt/ostream.h>
#include <iostream>
#include <string>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * formatted output using <a href="https://github.com/fmtlib/fmt">fmtlib</a>
 *
 * @param out
 * a @c FILE* pointer or an @c std::ostream object
 *
 * @param fmtstr
 * format string
 *
 * @param args
 * values to be printed out
 */
template <class Out, class Fmt, class... Ts>
void print (Out& out, const Fmt& fmtstr, const Ts&... args)
{
   fmt::print (out, fmtstr, args...);
}

/**
 * @brief
 * write the formatted output to an @c std::string object
 *
 * @param fmtstr
 * format string
 *
 * @param args
 * values to be printed out
 *
 * @return
 * the formatted @c std::string output
 */
template <class Fmt, class... Ts>
std::string format (const Fmt& fmtstr, const Ts&... args)
{
   return fmt::format (fmtstr, args...);
}
TINKER_NAMESPACE_END

#endif
