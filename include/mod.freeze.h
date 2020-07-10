#pragma once
#include "macro.h"
#include "mdprec.h"


namespace tinker {
TINKER_EXTERN real rateps;
TINKER_EXTERN int nrat;
TINKER_EXTERN int (*irat)[2];
TINKER_EXTERN real* krat;


TINKER_EXTERN int nratmol;
TINKER_EXTERN int (*iratmol)[2];


TINKER_EXTERN pos_prec* rattle_xold;
TINKER_EXTERN pos_prec* rattle_yold;
TINKER_EXTERN pos_prec* rattle_zold;
TINKER_EXTERN int* rattle_moved;
TINKER_EXTERN int* rattle_update;
TINKER_EXTERN int* rattle_notdone;
}
