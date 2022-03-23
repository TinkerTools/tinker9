#pragma once
#include "mod/rattle.h"
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

void hcKinetic();
void hcVirial();

void hcCenterOfMass(const pos_prec* atomx, const pos_prec* atomy, const pos_prec* atomz,
   pos_prec* molx, pos_prec* moly, pos_prec* molz);

// vi += s Vu
void hcVelIso(vel_prec s);
// vi += MatS Vu
void hcVelAn(vel_prec s[3][3]);
}
