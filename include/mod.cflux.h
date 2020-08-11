#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
TINKER_EXTERN real* bflx;
TINKER_EXTERN real* aflx;
TINKER_EXTERN real* abflx;


TINKER_EXTERN count_buffer naflx;
TINKER_EXTERN count_buffer nbflx;
TINKER_EXTERN energy_buffer pot;
}
