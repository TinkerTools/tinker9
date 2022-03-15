#pragma once
#include "macro.h"
#include <tinker/detail/ktrtor.hh>

namespace tinker {
// of size maxntt
TINKER_EXTERN int* tnx;
TINKER_EXTERN int* tny;
// of size (maxtgrd,maxntt) i.e. [maxntt][maxtgrd]
TINKER_EXTERN real (*ttx)[ktrtor::maxtgrd];
TINKER_EXTERN real (*tty)[ktrtor::maxtgrd];
// of size (maxtgrd2,maxntt) i.e. [maxntt][maxtgrd2]
TINKER_EXTERN real (*tbf)[ktrtor::maxtgrd2];
TINKER_EXTERN real (*tbx)[ktrtor::maxtgrd2];
TINKER_EXTERN real (*tby)[ktrtor::maxtgrd2];
TINKER_EXTERN real (*tbxy)[ktrtor::maxtgrd2];
}
