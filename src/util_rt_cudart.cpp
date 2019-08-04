#include "util_rt.h"

#ifndef TINKER_HOSTONLY
TINKER_NAMESPACE_BEGIN
void copyin_bytes(void* dst, const void* src, size_t count) {
  check_rt(cudaMemcpy(dst, src, count, cudaMemcpyHostToDevice));
}

void copyout_bytes(void* dst, const void* src, size_t count) {
  check_rt(cudaMemcpy(dst, src, count, cudaMemcpyDeviceToHost));
}

void copy_bytes(void* dst, const void* src, size_t count) {
  check_rt(cudaMemcpy(dst, src, count, cudaMemcpyDeviceToDevice));
}

void dealloc_stream(Stream s) { check_rt(cudaStreamDestroy(s)); }
void alloc_stream(Stream* s) { check_rt(cudaStreamCreate(s)); }
void sync_stream(Stream s) { check_rt(cudaStreamSynchronize(s)); }
void copy_bytes_async(void* dst, const void* src, size_t count, Stream s) {
  check_rt(cudaMemcpyAsync(dst, src, count, cudaMemcpyDeviceToDevice, s));
}
TINKER_NAMESPACE_END
#endif
