#include "ff/atom.h"
#include "ff/elec.h"
#include "ff/nblist.h"

namespace tinker {
extern void dfieldNonEwald_acc(real (*field)[3], real (*fieldp)[3]);
extern void dfieldNonEwald_cu(real (*field)[3], real (*fieldp)[3]);
void dfieldNonEwald(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfieldNonEwald_cu(field, fieldp);
   else
#endif
      dfieldNonEwald_acc(field, fieldp);
}
}

namespace tinker {
extern void dfieldEwaldRecipSelf_acc(real (*field)[3]);
static void dfieldEwaldRecipSelf(real (*field)[3], real (*fieldp)[3])
{
   dfieldEwaldRecipSelf_acc(field);
   darray::copy(g::q0, n, fieldp, field);
}

void dfieldEwaldReal_acc(real (*field)[3], real (*fieldp)[3]);
void dfieldEwaldReal_cu(real (*field)[3], real (*fieldp)[3]);
static void dfieldEwaldReal(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfieldEwaldReal_cu(field, fieldp);
   else
#endif
      dfieldEwaldReal_acc(field, fieldp);
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
extern void ufieldNonEwald_acc(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3]);
extern void ufieldNonEwald_cu(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3]);
void ufieldNonEwald(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufieldNonEwald_cu(uind, uinp, field, fieldp);
   else
#endif
      ufieldNonEwald_acc(uind, uinp, field, fieldp);
}
}

namespace tinker {
extern void ufieldEwaldRecipSelf_acc(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3]);
static void ufieldEwaldRecipSelf(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3])
{
   ufieldEwaldRecipSelf_acc(uind, uinp, field, fieldp);
}

extern void ufieldEwaldReal_acc(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3]);
extern void ufieldEwaldReal_cu(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3]);
void ufieldEwaldReal(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufieldEwaldReal_cu(uind, uinp, field, fieldp);
   else
#endif
      ufieldEwaldReal_acc(uind, uinp, field, fieldp);
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
