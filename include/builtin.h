#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
namespace builtin {
// int clz(int);
// int clz(long long);


namespace gcc {
inline int clz(int v)
{
   return __builtin_clz(v);
}


inline int clz(long long v)
{
   return __builtin_clzll(v);
}
}


using namespace gcc;
}
TINKER_NAMESPACE_END
