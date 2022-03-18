#pragma once
#include "precision.h"
#include "tool/rcman.h"

extern "C"
{
   struct RATTLE
   {};
   struct SHAKE
   {};
}

namespace tinker {
bool useRattle();
void rattleData(RcOp);
void rattle(time_prec dt, const pos_prec* xold, const pos_prec* yold, const pos_prec* zold);
void rattle2(time_prec dt, bool do_v);
void shake(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew, //
   const pos_prec* xold, const pos_prec* yold, const pos_prec* zold);
}

#include "glob.rattle.h"
