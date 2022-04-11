#include "ff/atom.h"
#include "ff/elec.h"
#include "ff/nblist.h"
#include "tool/externfunc.h"

namespace tinker {
TINKER_F2VOID(cu, 1, acc, 1, dfieldNonEwald, real (*)[3], real (*)[3]);
void dfieldNonEwald(real (*field)[3], real (*fieldp)[3])
{
   TINKER_F2CALL(cu, 1, acc, 1, dfieldNonEwald, field, fieldp);
}
}

namespace tinker {
TINKER_F2VOID(cu, 0, acc, 1, dfieldEwaldRecipSelf, real (*)[3]);
static void dfieldEwaldRecipSelf(real (*field)[3], real (*fieldp)[3])
{
   TINKER_F2CALL(cu, 0, acc, 1, dfieldEwaldRecipSelf, field);
   darray::copy(g::q0, n, fieldp, field);
}

TINKER_F2VOID(cu, 1, acc, 1, dfieldEwaldReal, real (*)[3], real (*)[3]);
static void dfieldEwaldReal(real (*field)[3], real (*fieldp)[3])
{
   TINKER_F2CALL(cu, 1, acc, 1, dfieldEwaldReal, field, fieldp);
}

void dfieldEwald(real (*field)[3], real (*fieldp)[3])
{
   dfieldEwaldRecipSelf(field, fieldp);
   dfieldEwaldReal(field, fieldp);
}
}

namespace tinker {
void dfield(real (*field)[3], real (*fieldp)[3])
{
   if (useEwald())
      dfieldEwald(field, fieldp);
   else
      dfieldNonEwald(field, fieldp);
}
}

namespace tinker {
TINKER_F2VOID(
   cu, 1, acc, 1, ufieldNonEwald, const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
void ufieldNonEwald(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   TINKER_F2CALL(cu, 1, acc, 1, ufieldNonEwald, uind, uinp, field, fieldp);
}
}

namespace tinker {
TINKER_F2VOID(cu, 0, acc, 1, ufieldEwaldRecipSelf, const real (*)[3], const real (*)[3],
   real (*)[3], real (*)[3]);
static void ufieldEwaldRecipSelf(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3])
{
   TINKER_F2CALL(cu, 0, acc, 1, ufieldEwaldRecipSelf, uind, uinp, field, fieldp);
}

TINKER_F2VOID(
   cu, 1, acc, 1, ufieldEwaldReal, const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
void ufieldEwaldReal(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3])
{
   TINKER_F2CALL(cu, 1, acc, 1, ufieldEwaldReal, uind, uinp, field, fieldp);
}

void ufieldEwald(const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   ufieldEwaldRecipSelf(uind, uinp, field, fieldp);
   ufieldEwaldReal(uind, uinp, field, fieldp);
}
}

namespace tinker {
void ufield(const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   if (useEwald())
      ufieldEwald(uind, uinp, field, fieldp);
   else
      ufieldNonEwald(uind, uinp, field, fieldp);
}
}
