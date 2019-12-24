#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
extern real* pinned_real64;
extern real* dptr_real64;


template <class A, class T1, class T2>
void dotprod(A* ans, const T1* a, const T2* b, int nelem);
TINKER_NAMESPACE_END
