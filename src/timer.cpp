#include "timer.h"
#include "io_print.h"
#include <cassert>
#include <chrono>
#include <vector>

TINKER_NAMESPACE_BEGIN
class Stopwatch
{
private:
   typedef std::chrono::steady_clock ClockType;
   std::vector<ClockType::time_point> s_;

public:
   void start();
   void lap(std::string log = "");
   void stop();
   void reset();
};

void Stopwatch::start()
{
   assert(s_.size() == 0);
   s_.emplace_back(ClockType::now());
}

void Stopwatch::lap(std::string log)
{
   auto now = ClockType::now();
   auto dur = now - s_.back();
   s_.emplace_back(ClockType::now());

   auto p = [=](double& v, std::string& u) {
      auto ns_i = std::chrono::nanoseconds(dur).count();
      double ns = static_cast<double>(ns_i);
      if (ns < 1000) {
         v = ns;
         u = "ns";
         return;
      }

      double us = ns / 1000;
      if (us < 1000) {
         v = us;
         u = "us";
         return;
      }

      double ms = us / 1000;
      if (ms < 1000) {
         v = ms;
         u = "ms";
         return;
      }

      double sec = ms / 1000;
      v = sec;
      u = " s";
      return;
   };

   double val;
   std::string u;
   p(val, u);
   if (log == "") {
      print(stdout, " {:>12.4f} {}           Lap {:<d}\n", val, u,
            s_.size() - 1);
   } else {
      print(stdout, " {:>12.4f} {}           {}\n", val, u, log);
   }
}

void Stopwatch::stop() {}

void Stopwatch::reset()
{
   s_.clear();
}

static Stopwatch sw_;

void stopwatch_start()
{
   sw_.start();
}

void stopwatch_lap(std::string log)
{
   sw_.lap(log);
}

void stopwatch_stop()
{
   sw_.stop();
}

void stopwatch_reset()
{
   sw_.reset();
}
TINKER_NAMESPACE_END
