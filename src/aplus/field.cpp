#include "ff/elec.h"
#include "ff/hippo/induce.h"
#include "ff/nblist.h"

namespace tinker {
extern void dfieldEwaldRecipSelf_acc(real (*field)[3]);
static void dfieldAplusEwaldRecipSelf(real (*field)[3])
{
   dfieldEwaldRecipSelf_acc(field);
}

extern void dfieldAplusEwaldReal_acc(real (*field)[3]);
extern void dfieldAplusEwaldReal_cu(real (*field)[3]);
static void dfieldAplusEwaldReal(real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfieldAplusEwaldReal_cu(field);
   else
#endif
      dfieldAplusEwaldReal_acc(field);
}

static void dfieldAplusEwald(real (*field)[3])
{
   dfieldAplusEwaldRecipSelf(field);
   dfieldAplusEwaldReal(field);
}
}

namespace tinker {
extern void dfieldAplusNonEwald_acc(real (*field)[3]);
extern void dfieldAplusNonEwald_cu(real (*field)[3]);
static void dfieldAplusNonEwald(real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfieldAplusNonEwald_cu(field);
   else
#endif
      dfieldAplusNonEwald_acc(field);
}

void dfieldAplus(real (*field)[3])
{
   if (useEwald())
      dfieldAplusEwald(field);
   else
      dfieldAplusNonEwald(field);
}
}

namespace tinker {
static void ufieldAplusEwaldRecipSelf(const real (*uind)[3], real (*field)[3])
{
   ufieldChgpenEwaldRecipSelf(uind, field);
}

extern void ufieldAplusEwaldReal_acc(const real (*uind)[3], real (*field)[3]);
extern void ufieldAplusEwaldReal_cu(const real (*uind)[3], real (*field)[3]);
static void ufieldAplusEwaldReal(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufieldAplusEwaldReal_cu(uind, field);
   else
#endif
      ufieldAplusEwaldReal_acc(uind, field);
}

static void ufieldAplusEwald(const real (*uind)[3], real (*field)[3])
{
   ufieldAplusEwaldRecipSelf(uind, field);
   ufieldAplusEwaldReal(uind, field);
}
}

namespace tinker {
extern void ufieldAplusNonEwald_acc(const real (*uind)[3], real (*field)[3]);
extern void ufieldAplusNonEwald_cu(const real (*uind)[3], real (*field)[3]);
static void ufieldAplusNonEwald(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufieldAplusNonEwald_cu(uind, field);
   else
#endif
      ufieldAplusNonEwald_acc(uind, field);
}

void ufieldAplus(const real (*uind)[3], real (*field)[3])
{
   if (useEwald())
      ufieldAplusEwald(uind, field);
   else
      ufieldAplusNonEwald(uind, field);
}
}
