#pragma once
#include "macro.h"
#include <cstring>
#include <type_traits>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup mem
 * \brief
 * A helper class to get the traits of the given 1D or 2D pointer type.
 *
 * Example:
 * | PTR        | deduce_ptr<PTR>::type | deduce_ptr<PTR>::n |
 * |:----------:|:------:|:---:|
 * | float*     | float  | 1   |
 * | int (*)[3] | int    | 3   |
 *
 * \tparam PTR
 * Must be a pointer type.
 */
template <class PTR>
struct deduce_ptr;


template <class T>
struct deduce_ptr<T*>
{
   typedef T type;
   static constexpr size_t n = 1;
};


template <class T, size_t N>
struct deduce_ptr<T (*)[N]>
{
   typedef T type;
   static constexpr size_t n = N;
};
TINKER_NAMESPACE_END
