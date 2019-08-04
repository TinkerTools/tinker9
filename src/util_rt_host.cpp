#include "util_rt.h"
#include <thread>

#ifdef TINKER_HOST
TINKER_NAMESPACE_BEGIN
void zero_bytes(void* ptr, size_t nbytes) { std::memset(ptr, 0, nbytes); }

void dealloc_bytes(void* ptr) { std::free(ptr); }

void alloc_bytes(void** ptr, size_t nbytes) { *ptr = std::malloc(nbytes); }

void copyin_bytes(void* dst, const void* src, size_t nbytes) {
  std::memcpy(dst, src, nbytes);
}

void copyout_bytes(void* dst, const void* src, size_t nbytes) {
  std::memcpy(dst, src, nbytes);
}

void copy_bytes(void* dst, const void* src, size_t nbytes) {
  std::memcpy(dst, src, nbytes);
}

void StreamSt::clear_front_() {
  auto exec = [&]() {
    auto& func = q_.front();
    func();
    mq_.lock();
    q_.pop();
    mq_.unlock();

    mi_.lock();
    idle_ = true;
    cvi_.notify_all();
    mi_.unlock();
  };
  std::thread(exec).detach();
}

void StreamSt::sync() {
  std::unique_lock<std::mutex> lck_idle(mi_);
  cvi_.wait(lck_idle, [=]() { return idle_; });
  while (q_.size()) {
    auto& func = q_.front();
    func();
    q_.pop();
  }
}

StreamSt::StreamSt()
    : mq_()
    , mi_()
    , cvi_()
    , q_()
    , idle_(true) {}

void dealloc_stream(Stream s) { delete s; }

void alloc_stream(Stream* ps) { *ps = new StreamSt; }

void sync_stream(Stream s) { s->sync(); }

void copy_bytes_async(void* dst, const void* src, size_t nbytes, Stream s) {
  s->add_async_call(std::memcpy, dst, src, nbytes);
}
TINKER_NAMESPACE_END
#endif
