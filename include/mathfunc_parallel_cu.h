#pragma once
#include "macro.h"
#include <cstring>


TINKER_NAMESPACE_BEGIN
namespace platform {
namespace cu {
template <class T>
void dotprod(T* ans, const T* a, const T* b, int nelem, int sync)
#if TINKER_HOST
{}
#endif
;
}
}
TINKER_NAMESPACE_END
