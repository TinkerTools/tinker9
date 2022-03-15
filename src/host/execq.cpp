#include "execq.h"
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>

namespace tinker {
class ExecQ::Impl
{};

void ExecQ::deallocate() {}

void ExecQ::allocate() {}

void ExecQ::begin_copyout() {}

void ExecQ::end_copyout() {}
}
