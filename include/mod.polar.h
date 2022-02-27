#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
TINKER_EXTERN real* polarity;
TINKER_EXTERN real* thole;
TINKER_EXTERN real* pdamp;
TINKER_EXTERN real* dirdamp;


TINKER_EXTERN real (*udir)[3];
TINKER_EXTERN real (*udirp)[3];
TINKER_EXTERN real (*uind)[3];
TINKER_EXTERN real (*uinp)[3];


//====================================================================//


TINKER_EXTERN count_buffer nep;
TINKER_EXTERN energy_buffer ep;
TINKER_EXTERN virial_buffer vir_ep;
TINKER_EXTERN grad_prec* depx;
TINKER_EXTERN grad_prec* depy;
TINKER_EXTERN grad_prec* depz;
TINKER_EXTERN energy_prec energy_ep;
TINKER_EXTERN virial_prec virial_ep[9];


TINKER_EXTERN real* polarity_inv;


TINKER_EXTERN real (*ufld)[3];
TINKER_EXTERN real (*dufld)[6];


TINKER_EXTERN real (*work01_)[3];
TINKER_EXTERN real (*work02_)[3];
TINKER_EXTERN real (*work03_)[3];
TINKER_EXTERN real (*work04_)[3];
TINKER_EXTERN real (*work05_)[3];
TINKER_EXTERN real (*work06_)[3];
TINKER_EXTERN real (*work07_)[3];
TINKER_EXTERN real (*work08_)[3];
TINKER_EXTERN real (*work09_)[3];
TINKER_EXTERN real (*work10_)[3];
}
