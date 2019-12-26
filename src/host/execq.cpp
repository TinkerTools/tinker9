#include "execq.h"
#include <condition_variable>
#include <cstdlib>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>


TINKER_NAMESPACE_BEGIN
class ExecQ::Impl
{
private:
   std::mutex mq_, mi_;
   std::condition_variable cvi_;
   std::queue<std::function<void()>> q_;
   bool idle_;


   void clear_front_()
   {
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


public:
   template <class F, class... Args>
   void add_async_call(F&& call, Args&&... args)
   {
      mq_.lock();
      q_.emplace(std::bind(call, args...));
      mq_.unlock();
      if (idle_) {
         idle_ = false;
         clear_front_();
      }
   }


   Impl()
      : mq_()
      , mi_()
      , cvi_()
      , q_()
      , idle_(true)
   {}


   void sync()
   {
      std::unique_lock<std::mutex> lck_idle(mi_);
      cvi_.wait(lck_idle, [=]() { return idle_; });
      while (q_.size()) {
         auto& func = q_.front();
         func();
         q_.pop();
      }
   }
};


void ExecQ::deallocate()
{
   delete ptr;
}


void ExecQ::allocate()
{
   ptr = new ExecQ::Impl;
}


void ExecQ::synchronize()
{
   ptr->sync();
}


void ExecQ::copy_bytes(void* dst, const void* src, size_t nbytes)
{
   ptr->add_async_call(std::memcpy, dst, src, nbytes);
}
TINKER_NAMESPACE_END
