#ifndef TINKER_UTIL_RC_MAN_H_
#define TINKER_UTIL_RC_MAN_H_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN
void fortran_runtime_initialize(int, char**);
void fortran_runtime_finish();
TINKER_NAMESPACE_END

#endif
