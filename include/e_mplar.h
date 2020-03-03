#pragma once
#include "dev_array.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int nmdpuexclude;
TINKER_EXTERN pointer<int, 2> mdpuexclude;
TINKER_EXTERN pointer<real, 4> mdpuexclude_scale;


void emplar_data(rc_op op);
void emplar_cu(int vers);
TINKER_NAMESPACE_END
