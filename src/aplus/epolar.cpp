#include "ff/amoeba/epolar.h"
#include "ff/amoebamod.h"
#include "ff/hippo/epolar.h"
#include "ff/hippo/induce.h"
#include "ff/nblist.h"
#include "tool/externfunc.h"

namespace tinker {
TINKER_F2VOID(cu, 1, acc, 1, epolarAplusEwaldReal, int, int, const real (*)[3]);
static void epolarAplusEwaldReal(int vers, int use_cf)
{
   TINKER_F2CALL(cu, 1, acc, 1, epolarAplusEwaldReal, vers, use_cf, uind);
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

   induce2(uind);
   if (edot)
      epolar0DotProd(uind, udir);
   if (vers != calc::v0) {
      epolarAplusEwaldReal(ver2, use_cf);
      epolarAplusEwaldRecipSelf(ver2, use_cf);
   }
}
}

namespace tinker {
TINKER_F2VOID(cu, 1, acc, 1, epolarAplusNonEwald, int, int, const real (*)[3]);
void epolarAplusNonEwald(int vers, int use_cf)
{
   bool edot = vers & calc::energy; // if not do_e, edot = false
   if (vers & calc::energy and vers & calc::analyz)
      edot = false; // if do_e and do_a, edot = false
   int ver2 = vers;
   if (edot)
      ver2 &= ~calc::energy; // toggle off the calc::energy flag

   induce2(uind);
   if (edot)
      epolar0DotProd(uind, udir);
   if (vers != calc::v0)
      TINKER_F2CALL(cu, 1, acc, 1, epolarAplusNonEwald, ver2, use_cf, uind);
}
}
