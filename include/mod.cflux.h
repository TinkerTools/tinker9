#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
TINKER_EXTERN real* bflx;
TINKER_EXTERN real (*aflx)[2];
TINKER_EXTERN real (*abflx)[2];
TINKER_EXTERN int* atomic;
TINKER_EXTERN int (*balist)[2];

//TINKER_EXTERN count_buffer naflx;
//TINKER_EXTERN count_buffer nbflx;

TINKER_EXTERN real* pdelta;
TINKER_EXTERN real* pot;
TINKER_EXTERN real* decfx;
TINKER_EXTERN real* decfy;
TINKER_EXTERN real* decfz;

TINKER_EXTERN int (*couple_i12)[8];
TINKER_EXTERN int* couple_n12;

}
