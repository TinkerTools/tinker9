#include "ff/hippo/empole.h"
#include "ff/nblist.h"

namespace tinker {
extern void empoleAplusEwaldRealSelf_acc(int vers, int useCF);
extern void empoleAplusEwaldRealSelf_cu(int vers, int useCF);
static void empoleAplusEwaldRealSelf(int vers, int useCF)
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      empoleAplusEwaldRealSelf_cu(vers, useCF);
   else
#endif
      empoleAplusEwaldRealSelf_acc(vers, useCF);
}

static void empoleAplusEwaldRecip(int vers, int useCF)
{
   empoleChgpenEwaldRecip(vers, useCF);
}

void empoleAplusEwald(int vers, int useCF)
{
   empoleAplusEwaldRealSelf(vers, useCF);
   empoleAplusEwaldRecip(vers, useCF);
}
}

namespace tinker {
extern void empoleAplusNonEwald_acc(int, int);
extern void empoleAplusNonEwald_cu(int, int);
void empoleAplusNonEwald(int vers, int useCF)
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      empoleAplusNonEwald_cu(vers, useCF);
   else
#endif
      empoleAplusNonEwald_acc(vers, useCF);
}
}
