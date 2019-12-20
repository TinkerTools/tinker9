#pragma once
#include "dev_array.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int nmdpuexclude;
TINKER_EXTERN device_pointer<int, 2> mdpuexclude;
TINKER_EXTERN device_pointer<real, 4> mdpuexclude_scale;


void emplar_cu(int vers);
void emplar_data(rc_op op);
TINKER_NAMESPACE_END
