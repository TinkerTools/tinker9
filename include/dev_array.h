#ifndef TINKER_DEV_ARRAY_H_
#define TINKER_DEV_ARRAY_H_

#include "dev_memory.h"

TINKER_NAMESPACE_BEGIN
/**
 * \typedef device_array
 * \ingroup util
 * A collection of device array related functions.
 * \see DeviceArray
 * \see DeviceMemory
 */
typedef DeviceArray<DeviceMemory> device_array;

/**
 * \typedef device_pointer
 * \ingroup util
 * Based on the template parameters \c T and \c N, this type is either
 * defined to `T*` or `T(*)[N]` when \c N is greater than 1.
 * \c N is set to 1 by default.
 */
template <class T, size_t N = 1>
using device_pointer = typename device_array::ptr<T, N>::type;
TINKER_NAMESPACE_END

#endif
