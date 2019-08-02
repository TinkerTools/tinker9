#ifndef TINKER_HOSTONLY

#  include "util_rt.h"

TINKER_NAMESPACE_BEGIN
void dealloc_stream(Stream s) { check_rt(cudaStreamDestroy(s)); }

void alloc_stream(Stream* s) { check_rt(cudaStreamCreate(s)); }

void sync_stream(Stream s) { check_rt(cudaStreamSynchronize(s)); }

void copy_memory(void* dst, const void* src, size_t count, CopyDirection k) {
  auto kind = cudaMemcpyHostToHost;
  if (k == CopyDirection::HostToDevice)
    kind = cudaMemcpyHostToDevice;
  else if (k == CopyDirection::DeviceToHost)
    kind = cudaMemcpyDeviceToHost;
  else if (k == CopyDirection::DeviceToDevice)
    kind = cudaMemcpyDeviceToDevice;
  else
    assert(false);

  check_rt(cudaMemcpy(dst, src, count, kind));
}

void copy_memory_async(void* dst, const void* src, size_t count,
                       CopyDirection k, Stream s) {
  auto kind = cudaMemcpyHostToHost;
  if (k == CopyDirection::HostToDevice)
    kind = cudaMemcpyHostToDevice;
  else if (k == CopyDirection::DeviceToHost)
    kind = cudaMemcpyDeviceToHost;
  else if (k == CopyDirection::DeviceToDevice)
    kind = cudaMemcpyDeviceToDevice;
  else
    assert(false);

  check_rt(cudaMemcpyAsync(dst, src, count, kind, s));
}

TINKER_NAMESPACE_END

#endif
