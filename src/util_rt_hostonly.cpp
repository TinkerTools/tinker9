#include "util_rt.h"
#include <thread>

#ifdef TINKER_HOSTONLY
TINKER_NAMESPACE_BEGIN
void copyin_bytes(void* dst, const void* src, size_t count) {
  std::memcpy(dst, src, count);
}

void copyout_bytes(void* dst, const void* src, size_t count) {
  std::memcpy(dst, src, count);
}

void copy_bytes(void* dst, const void* src, size_t count) {
  std::memcpy(dst, src, count);
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
void copy_bytes_async(void* dst, const void* src, size_t count, Stream s) {
  s->add_async_call(std::memcpy, dst, src, count);
}
TINKER_NAMESPACE_END
#endif
