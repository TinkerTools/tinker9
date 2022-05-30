#pragma once
#include <cstddef>
#include <type_traits>

namespace tinker {
/**
 * \ingroup cpp_syntax
 * A helper class to get the traits of the given 1D or 2D pointer type.
 *
 * Example:
 * | PTR        | PtrTrait<PTR>::type | PtrTrait<PTR>::n |
 * |:----------:|:-------------------:|:----------------:|
 * | float*     | float               | 1                |
 * | int (*)[3] | int                 | 3                |
 *
 * \tparam PTR  A pointer type.
 */
template <class PTR>
struct PtrTrait;

template <class T>
struct PtrTrait<T*>
{
   typedef T type;
   static constexpr size_t n = 1;
};

template <class T, size_t N>
struct PtrTrait<T (*)[N]>
{
   typedef T type;
   static constexpr size_t n = N;
};
}
