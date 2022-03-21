#pragma once
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
   void beginCopyout();
   void endCopyout();
};
}
