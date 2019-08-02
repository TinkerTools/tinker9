#ifndef TINKER_HOSTONLY

#  include "util_array.h"
#  include "util_rt.h"

TINKER_NAMESPACE_BEGIN
void dealloc_array(void* ptr) { check_rt(cudaFree(ptr)); }

void alloc_array(void** ptr, size_t nbytes) {
  check_rt(cudaMalloc(ptr, nbytes));
}

template <class T>
void zero_array_tmpl(T* dst, int nelem) {
  size_t size = sizeof(T) * nelem;
  check_rt(cudaMemset(dst, 0, size));
}

void zero_array(int* dst, int nelem) { zero_array_tmpl(dst, nelem); }

void zero_array(float* dst, int nelem) { zero_array_tmpl(dst, nelem); }

void zero_array(double* dst, int nelem) { zero_array_tmpl(dst, nelem); }
TINKER_NAMESPACE_END

#endif
