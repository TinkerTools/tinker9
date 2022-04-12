#include "ff/elec.h"
#include "ff/nblist.h"
#include "tool/externfunc.h"

namespace tinker {
TINKER_FVOID2(cu, 0, acc, 1, dfieldEwaldRecipSelf, real (*)[3]);
static void dfieldChgpenEwaldRecipSelf(real (*field)[3])
{
   TINKER_FCALL2(cu, 0, acc, 1, dfieldEwaldRecipSelf, field);
}

TINKER_FVOID2(cu, 1, acc, 1, dfieldChgpenEwaldReal, real (*)[3]);
static void dfieldChgpenEwaldReal(real (*field)[3])
{
   TINKER_FCALL2(cu, 1, acc, 1, dfieldChgpenEwaldReal, field);
}

static void dfieldChgpenEwald(real (*field)[3])
{
   dfieldChgpenEwaldRecipSelf(field);
   dfieldChgpenEwaldReal(field);
}
}

namespace tinker {
TINKER_FVOID2(cu, 1, acc, 1, dfieldChgpenNonEwald, real (*)[3]);
static void dfieldChgpenNonEwald(real (*field)[3])
{
   TINKER_FCALL2(cu, 1, acc, 1, dfieldChgpenNonEwald, field);
}

void dfieldChgpen(real (*field)[3])
{
   if (useEwald())
      dfieldChgpenEwald(field);
   else
      dfieldChgpenNonEwald(field);
}
}

namespace tinker {
TINKER_FVOID2(cu, 0, acc, 1, ufieldChgpenEwaldRecipSelf, const real (*)[3], real (*)[3]);
void ufieldChgpenEwaldRecipSelf(const real (*uind)[3], real (*field)[3])
{
   TINKER_FCALL2(cu, 0, acc, 1, ufieldChgpenEwaldRecipSelf, uind, field);
}

TINKER_FVOID2(cu, 1, acc, 1, ufieldChgpenEwaldReal, const real (*)[3], real (*)[3]);
static void ufieldChgpenEwaldReal(const real (*uind)[3], real (*field)[3])
{
   TINKER_FCALL2(cu, 1, acc, 1, ufieldChgpenEwaldReal, uind, field);
}

static void ufieldChgpenEwald(const real (*uind)[3], real (*field)[3])
{
   ufieldChgpenEwaldRecipSelf(uind, field);
   ufieldChgpenEwaldReal(uind, field);
}
}

namespace tinker {
TINKER_FVOID2(cu, 1, acc, 1, ufieldChgpenNonEwald, const real (*)[3], real (*)[3]);
static void ufieldChgpenNonEwald(const real (*uind)[3], real (*field)[3])
{
   TINKER_FCALL2(cu, 1, acc, 1, ufieldChgpenNonEwald, uind, field);
}

void ufieldChgpen(const real (*uind)[3], real (*field)[3])
{
   if (useEwald())
      ufieldChgpenEwald(uind, field);
   else
      ufieldChgpenNonEwald(uind, field);
}
}
