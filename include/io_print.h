#pragma once
#include "macro.h"
#include <fmt/ostream.h>
#include <iostream>
#include <string>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup io
 * \brief
 * Formatted output using <a href="https://github.com/fmtlib/fmt">fmtlib</a>.
 * \param out
 * A `FILE*` pointer or an `std::ostream` object.
 * \param fmtstr
 * Format string.
 * \param args
 * Values to be printed out.
 */
template <class Out, class Fmt, class... Ts>
void print(Out& out, const Fmt& fmtstr, const Ts&... args)
{
   fmt::print(out, fmtstr, args...);
}

/**
 * \ingroup io
 * \brief
 * Write the formatted output to an `std::string` object.
 * \param fmtstr
 * Format string.
 * \param args
 * Values to be printed out.
 * \return
 * The formatted `std::string` output.
 */
template <class Fmt, class... Ts>
std::string format(const Fmt& fmtstr, const Ts&... args)
{
   return fmt::format(fmtstr, args...);
}
TINKER_NAMESPACE_END
