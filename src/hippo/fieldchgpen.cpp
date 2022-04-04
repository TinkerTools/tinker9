#include "ff/elec.h"
#include "ff/nblist.h"

namespace tinker {
extern void dfieldEwaldRecipSelf_acc(real (*field)[3]);
static void dfieldChgpenEwaldRecipSelf(real (*field)[3])
{
   dfieldEwaldRecipSelf_acc(field);
}

extern void dfieldChgpenEwaldReal_acc(real (*field)[3]);
extern void dfieldChgpenEwaldReal_cu(real (*field)[3]);
static void dfieldChgpenEwaldReal(real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfieldChgpenEwaldReal_cu(field);
   else
#endif
      dfieldChgpenEwaldReal_acc(field);
}

static void dfieldChgpenEwald(real (*field)[3])
{
   dfieldChgpenEwaldRecipSelf(field);
   dfieldChgpenEwaldReal(field);
}
}

namespace tinker {
extern void dfieldChgpenNonEwald_acc(real (*field)[3]);
extern void dfieldChgpenNonEwald_cu(real (*field)[3]);
static void dfieldChgpenNonEwald(real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfieldChgpenNonEwald_cu(field);
   else
#endif
      dfieldChgpenNonEwald_acc(field);
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
extern void ufieldChgpenEwaldRecipSelf_acc(const real (*uind)[3], real (*field)[3]);
static void ufieldChgpenEwaldRecipSelf(const real (*uind)[3], real (*field)[3])
{
   ufieldChgpenEwaldRecipSelf_acc(uind, field);
}

extern void ufieldChgpenEwaldReal_acc(const real (*uind)[3], real (*field)[3]);
extern void ufieldChgpenEwaldReal_cu(const real (*uind)[3], real (*field)[3]);
static void ufieldChgpenEwaldReal(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufieldChgpenEwaldReal_cu(uind, field);
   else
#endif
      ufieldChgpenEwaldReal_acc(uind, field);
}

static void ufieldChgpenEwald(const real (*uind)[3], real (*field)[3])
{
   ufieldChgpenEwaldRecipSelf(uind, field);
   ufieldChgpenEwaldReal(uind, field);
}
}

namespace tinker {
extern void ufieldChgpenNonEwald_acc(const real (*uind)[3], real (*field)[3]);
extern void ufieldChgpenNonEwald_cu(const real (*uind)[3], real (*field)[3]);
static void ufieldChgpenNonEwald(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufieldChgpenNonEwald_cu(uind, field);
   else
#endif
      ufieldChgpenNonEwald_acc(uind, field);
}

void ufieldChgpen(const real (*uind)[3], real (*field)[3])
{
   if (useEwald())
      ufieldChgpenEwald(uind, field);
   else
      ufieldChgpenNonEwald(uind, field);
}
}
