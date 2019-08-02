#ifdef TINKER_HOSTONLY

#  include "util_array.h"

TINKER_NAMESPACE_BEGIN
void dealloc_array(void* ptr) { std::free(ptr); }

void alloc_array(void** ptr, size_t nbytes) { *ptr = std::malloc(nbytes); }

template <class T>
void zero_array_tmpl(T* dst, int nelem) {
  size_t size = sizeof(T) * nelem;
  std::memset(dst, 0, size);
}

void zero_array(int* dst, int nelem) { zero_array_tmpl(dst, nelem); }

void zero_array(float* dst, int nelem) { zero_array_tmpl(dst, nelem); }

void zero_array(double* dst, int nelem) { zero_array_tmpl(dst, nelem); }
TINKER_NAMESPACE_END

#endif
