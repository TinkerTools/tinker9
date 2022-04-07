#include "ff/amoeba/epolar.h"
#include "ff/amoebamod.h"
#include "ff/aplus/induce.h"
#include "ff/hippo/epolar.h"
#include "ff/nblist.h"

namespace tinker {
extern void epolarAplusEwaldReal_acc(int vers, int use_cf, const real (*d)[3]);
extern void epolarAplusEwaldReal_cu(int vers, int use_cf, const real (*d)[3]);
static void epolarAplusEwaldReal(int vers, int use_cf)
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      epolarAplusEwaldReal_cu(vers, use_cf, uind);
   else
#endif
      epolarAplusEwaldReal_acc(vers, use_cf, uind);
}

static void epolarAplusEwaldRecipSelf(int vers, int use_cf)
{
   epolarChgpenEwaldRecipSelf(vers, use_cf);
}

void epolarAplusEwald(int vers, int use_cf)
{
   bool edot = vers & calc::energy; // if not do_e, edot = false
   if (vers & calc::energy and vers & calc::analyz)
      edot = false; // if do_e and do_a, edot = false
   int ver2 = vers;
   if (edot)
      ver2 &= ~calc::energy; // toggle off the calc::energy flag

   induce3(uind);
   if (edot)
      epolar0DotProd(uind, udir);
   if (vers != calc::v0) {
      epolarAplusEwaldReal(ver2, use_cf);
      epolarAplusEwaldRecipSelf(ver2, use_cf);
   }
}
}

namespace tinker {
extern void epolarAplusNonEwald_acc(int vers, int use_cf, const real (*d)[3]);
extern void epolarAplusNonEwald_cu(int vers, int use_cf, const real (*d)[3]);
void epolarAplusNonEwald(int vers, int use_cf)
{
   bool edot = vers & calc::energy; // if not do_e, edot = false
   if (vers & calc::energy and vers & calc::analyz)
      edot = false; // if do_e and do_a, edot = false
   int ver2 = vers;
   if (edot)
      ver2 &= ~calc::energy; // toggle off the calc::energy flag

   induce3(uind);
   if (edot)
      epolar0DotProd(uind, udir);
   if (vers != calc::v0) {
#if TINKER_CUDART
      if (mlistVersion() & Nbl::SPATIAL)
         epolarAplusNonEwald_cu(ver2, use_cf, uind);
      else
#endif
         epolarAplusNonEwald_acc(ver2, use_cf, uind);
   }
}
}
