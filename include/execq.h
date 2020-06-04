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
   void synchronize();
   void copy_bytes(void* dst, const void* src, size_t nbytes);
};
}
