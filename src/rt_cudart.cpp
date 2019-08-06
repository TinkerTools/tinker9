#include "rt.h"
#include "error.h"

#ifndef TINKER_HOST
TINKER_NAMESPACE_BEGIN
void zero_bytes(void* dst, size_t nbytes) {
  check_rt(cudaMemset(dst, 0, nbytes));
}

void dealloc_bytes(void* ptr) { check_rt(cudaFree(ptr)); }

void alloc_bytes(void** ptr, size_t nbytes) {
  check_rt(cudaMalloc(ptr, nbytes));
}

void copyin_bytes(void* dst, const void* src, size_t nbytes) {
  check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyHostToDevice));
}

void copyout_bytes(void* dst, const void* src, size_t nbytes) {
  check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToHost));
}

void copy_bytes(void* dst, const void* src, size_t nbytes) {
  check_rt(cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToDevice));
}

void dealloc_stream(Stream s) { check_rt(cudaStreamDestroy(s)); }

void alloc_stream(Stream* s) { check_rt(cudaStreamCreate(s)); }

void sync_stream(Stream s) { check_rt(cudaStreamSynchronize(s)); }

void copy_bytes_async(void* dst, const void* src, size_t nbytes, Stream s) {
  check_rt(cudaMemcpyAsync(dst, src, nbytes, cudaMemcpyDeviceToDevice, s));
}
TINKER_NAMESPACE_END
#endif
