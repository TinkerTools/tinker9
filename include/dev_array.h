#ifndef TINKER_DEV_ARRAY_H_
#define TINKER_DEV_ARRAY_H_

#include "dev_memory.h"

TINKER_NAMESPACE_BEGIN
typedef DeviceArray<DeviceMemory> device_array;

template <class T, size_t N = 1>
using device_pointer = typename device_array::ptr<T, N>::type;
TINKER_NAMESPACE_END

#endif
