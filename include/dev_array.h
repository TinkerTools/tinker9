#pragma once
#include "dev_memory.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup mem
 * A collection of device array related functions.
 * \see DeviceArray
 * \see DeviceMemory
 */
typedef DeviceArray<DeviceMemory> device_array;


/**
 * \ingroup mem
 * Based on the template parameters `T` and `N`, this type is either
 * defined to `T*` or `T(*)[N]` when `N` is greater than 1.
 * `N` is set to 1 by default.
 */
template <class T, size_t N = 1>
using device_pointer = typename device_array::ptr<T, N>::type;
TINKER_NAMESPACE_END
