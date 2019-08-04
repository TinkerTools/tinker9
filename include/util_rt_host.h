#ifndef TINKER_UTIL_RT_HOST_H_
#define TINKER_UTIL_RT_HOST_H_

#include "macro.h"
#include <fftw3.h>

#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>

TINKER_NAMESPACE_BEGIN
class StreamSt {
private:
  std::mutex mq_, mi_;
  std::condition_variable cvi_;
  std::queue<std::function<void()>> q_;
  bool idle_;
  void clear_front_();

public:
  template <class F, class... Args>
  void add_async_call(F&& call, Args&&... args) {
    mq_.lock();
    q_.emplace(std::bind(call, args...));
    mq_.unlock();
    if (idle_) {
      idle_ = false;
      clear_front_();
    }
  }
  void sync();
  StreamSt();
};
typedef StreamSt* Stream;

struct FFTPlan {
#if defined(TINKER_SINGLE_PRECISION)
  fftwf_plan
#elif defined(TINKER_DOUBLE_PRECISION)
  fftw_plan
#else
  static_assert(false, "");
#endif
      planf, ///< fft front plan
      planb; ///< fft back plan
};
TINKER_NAMESPACE_END

#endif
