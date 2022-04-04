#include "ff/amoeba/elecamoeba.h"
#include "ff/amoeba/empole.h"
#include "ff/energy.h"
#include "math/zero.h"
#include "tool/error.h"

namespace tinker {
extern void emplar_cu(int);
void emplar(int vers)
{
#if TINKER_CUDART
   bool do_v = vers & calc::virial;

   zeroOnHost(energy_em, virial_em);

   mpoleInit(vers);
   emplar_cu(vers);
   torque(vers, demx, demy, demz);
   if (do_v) {
      VirialBuffer u2 = vir_trq;
      virial_prec v2[9];
      virialReduce(v2, u2);
      for (int iv = 0; iv < 9; ++iv)
         virial_elec[iv] += v2[iv];
   }
#else
   (void)vers;
   TINKER_THROW("EMPLAR is only available for PBC systems in CUDA.\n");
#endif
}
}
