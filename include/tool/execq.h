#pragma once

namespace tinker {
inline namespace v1 {
/// \ingroup async
/// \brief A helper class for mdsave.
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
}
