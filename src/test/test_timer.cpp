#include "test.h"
#include "timer.h"
#include <thread>

using namespace TINKER_NAMESPACE;

TEST_CASE("Stopwatch", "[noassert]")
{
   stopwatch_start();

   std::this_thread::sleep_for(std::chrono::milliseconds(10));
   stopwatch_lap();

   std::this_thread::sleep_for(std::chrono::microseconds(150));
   stopwatch_lap("Lap 2 Info");

   std::this_thread::sleep_for(std::chrono::milliseconds(20));
   stopwatch_lap("Lap 3 Info");

   stopwatch_stop();
   stopwatch_reset();
}
