#pragma once
#include "macro.h"
#include "mdprec.h"


namespace tinker {
TINKER_EXTERN pos_prec rateps;


TINKER_EXTERN int nratwt;            // rattle water
TINKER_EXTERN int (*iratwt)[3];      // atoms a b c
TINKER_EXTERN pos_prec (*kratwt)[3]; // lengths ab ac bc


TINKER_EXTERN int nratch;       // rattle methine group
TINKER_EXTERN int (*iratch)[2]; // atoms a b (-C H)
TINKER_EXTERN pos_prec* kratch; // length ab


TINKER_EXTERN int nratch2;            // rattle methylene group
TINKER_EXTERN int (*iratch2)[3];      // atoms a b c (-C H H)
TINKER_EXTERN pos_prec (*kratch2)[2]; // lengths ab ac


TINKER_EXTERN int nratch3;            // rattle methyl group
TINKER_EXTERN int (*iratch3)[4];      // atoms a b c d (-C H H H)
TINKER_EXTERN pos_prec (*kratch3)[3]; // lengths ab ac ad


TINKER_EXTERN int nrat;
TINKER_EXTERN int (*irat)[2];
TINKER_EXTERN pos_prec* krat;


TINKER_EXTERN int nratmol;
TINKER_EXTERN int (*iratmol)[2];


TINKER_EXTERN pos_prec* rattle_xold;
TINKER_EXTERN pos_prec* rattle_yold;
TINKER_EXTERN pos_prec* rattle_zold;
}
