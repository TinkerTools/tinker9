# -*- coding: utf-8 -*-

'''
This script detects the compute capabilities (X.Y) of the the GPU cards and
returns them in a comma-separated XY string, if multiple GPUs are detected.

The auto-detected compute capabilities (XY) are used by the CMakeLists.txt only
if they are not explicitly set by the user. For the following situations:
    - the compiler does not support the XY value, e.g., XY=52;
    - the compiler is older than the device;
    - you don't have python3 in your system;
    - etc.
you must set the compute capabilities explicitly when cmake is configured.

Reference:
https://gist.github.com/f0k/63a664160d016a491b2cbea15913d549
by Jan Schlüter
'''


import ctypes
import sys


def CommaSeparatedCCString():
    unknown_cc = 'UnknownCC'

    libnames = ('libcuda.so', 'libcuda.dylib', 'nvcuda.dll', 'cuda.dll')
    for libname in libnames:
        try:
            cuda = ctypes.CDLL(libname)
        except OSError:
            continue
        else:
            break
    else:
        return unknown_cc

    ############################################################

    def checkCall(f, *args):
        CUDA_SUCCESS = 0
        result = ctypes.c_int()

        result = f(*args)
        if result != CUDA_SUCCESS:
            error_str = ctypes.c_char_p()
            cuda.cuGetErrorString(result, ctypes.byref(error_str))
            print('{} failed with error code {}: {}'.format(f.__name__, result, error_str.value.decode()))
            sys.exit(1)

    ############################################################

    cc_set = set()

    # from cuda.h
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76

    checkCall(cuda.cuInit, 0)

    nGpus = ctypes.c_int()
    checkCall(cuda.cuDeviceGetCount, ctypes.byref(nGpus))

    for i in range(nGpus.value):
        device = ctypes.c_int()
        checkCall(cuda.cuDeviceGet, ctypes.byref(device), i)

        cc_major = ctypes.c_int()
        cc_minor = ctypes.c_int()
        checkCall(cuda.cuDeviceGetAttribute, ctypes.byref(cc_major), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device)
        checkCall(cuda.cuDeviceGetAttribute, ctypes.byref(cc_minor), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device)
        cc = cc_major.value * 10 + cc_minor.value
        cc_set.add(cc)

    cc_list = []
    for c in sorted(cc_set):
        cc_list.append(c)

    cc_len = len(cc_list)
    if cc_len:
        cc_str = '{}'.format(cc_list[0])
        for i in range(1, cc_len):
            cc_str = cc_str + ',{}'.format(cc_list[i])
        return cc_str
    else:
        return unknown_cc


if __name__ == '__main__':
    print(CommaSeparatedCCString(), end='')
