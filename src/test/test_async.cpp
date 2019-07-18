#include "test/test.h"
#include <chrono>
#include <condition_variable>
#include <future>
#include <mutex>
#include <random>

using namespace std;
using namespace chrono;

static int gpu_value;
static mutex mtx_copy, mtx_write;
static condition_variable cv_copy, cv_write;
static bool idle_copy, idle_write;
static future<void> fut_copy_then_write;

static mutex mtx_cout;
static int write_ms, copy_ms;

static void init() {
  gpu_value = 0;
  idle_copy = false;
  idle_write = true;

  random_device rd;
  default_random_engine e(rd());
  uniform_int_distribution<int> ud(51, 300);
  copy_ms = ud(e);
  write_ms = ud(e);

  mtx_cout.lock();
  printf("RANDOM LAPSE: COPY %4d us; WRITE %4d us.\n\n", copy_ms, write_ms);
  mtx_cout.unlock();
}

static void write_step(int idx) {
  // time spent on copying
  this_thread::sleep_for(microseconds(copy_ms));
  int copyout_value = gpu_value;
  int host_value = copyout_value;
  mtx_copy.lock();
  idle_copy = true;
  cv_copy.notify_all();
  mtx_copy.unlock();

  mtx_cout.lock();
  printf("== ASYNC_WRITE_1 %4d | GPU VALUE %4d | HOST %4d\n", idx, gpu_value,
         host_value);
  mtx_cout.unlock();

  // time spent on writing
  this_thread::sleep_for(microseconds(write_ms));
  mtx_write.lock();
  idle_write = true;
  cv_write.notify_all();
  mtx_write.unlock();

  mtx_cout.lock();
  printf("== ASYNC_WRITE_2 %4d | GPU VALUE %4d | HOST %4d\n", idx, gpu_value,
         host_value);
  mtx_cout.unlock();
}

static void save_step(int i) {
  mtx_cout.lock();
  printf("SAVE %4d | %s\n", i, (idle_write ? "WRITE_IDLE" : "WRITE_BUSY"));
  mtx_cout.unlock();

  // in copy_then_write():
  // 1. async copy gpu value to the host;
  // 2. async write host value to external files

  // this routine must wait for the previous writing operations to finish before
  // creating a new async thread; in copy_then_write() idle_write is set to TRUE
  // and is then notified to other threads once writing is finished
  unique_lock<mutex> lck_write(mtx_write);
  cv_write.wait(lck_write, [=]() { return idle_write; });
  idle_write = false;

  mtx_cout.lock();
  printf("## BEFORE ASYNC_WRITE %3d\n", i);
  mtx_cout.unlock();

  fut_copy_then_write = async(launch::async, write_step, i);

  mtx_cout.lock();
  printf("## AFTER  ASYNC_WRITE %3d\n************\n\n", i);
  mtx_cout.unlock();

  // gpu_value may change after exiting this routine even BEFORE
  // copy_then_write() being called in the async thread; in copy_then_write()
  // idle_copy is set to TRUE and is then notified to other threads once copy is
  // finished
  unique_lock<mutex> lck_copy(mtx_copy);
  cv_copy.wait(lck_copy, [=]() { return idle_copy; });
  idle_copy = false;
}

static void sync_write() { fut_copy_then_write.get(); }

static void func() {
  init();

  const int maxstep = 12;
  const int iwrite = 3;
  for (int i = 1; i <= maxstep; ++i) {
    gpu_value -= 1;
    if (i % iwrite == 0)
      save_step(i);

    mtx_cout.lock();
    printf("STEP %4d | GPU VALUE %4d\n", i, gpu_value);
    mtx_cout.unlock();
  }
  sync_write();
}

TEST_CASE("Host-Async", "[noassert][host][async]") { func(); }
