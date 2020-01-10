#pragma once
#include "rc_man.h"
#include <cstddef>


TINKER_NAMESPACE_BEGIN
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
TINKER_NAMESPACE_END
