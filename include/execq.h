#pragma once
#include "tool/rc_man.h"
#include <cstddef>


namespace tinker {
class ExecQ
{
private:
   class Impl;
   Impl* ptr;


public:
   void deallocate();
   void allocate();
   void begin_copyout();
   void end_copyout();
};
}
