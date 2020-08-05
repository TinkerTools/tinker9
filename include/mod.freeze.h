#pragma once
#include "macro.h"
#include "mdprec.h"


namespace tinker {
TINKER_EXTERN pos_prec rateps;


TINKER_EXTERN int nratwt;            // rattle water
TINKER_EXTERN int (*iratwt)[3];      // atoms a b c
TINKER_EXTERN pos_prec (*kratwt)[3]; // lengths ab ac bc


TINKER_EXTERN int nrat;
TINKER_EXTERN int (*irat)[2];
TINKER_EXTERN pos_prec* krat;


TINKER_EXTERN int nratmol;
TINKER_EXTERN int (*iratmol)[2];


TINKER_EXTERN pos_prec* rattle_xold;
TINKER_EXTERN pos_prec* rattle_yold;
TINKER_EXTERN pos_prec* rattle_zold;


// These three variables are only used for the old rattle subroutines.
TINKER_EXTERN int* rattle_moved;
TINKER_EXTERN int* rattle_update;
TINKER_EXTERN int* rattle_bigdelta;
}
