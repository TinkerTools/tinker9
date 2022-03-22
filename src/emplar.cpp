#include "ff/amoeba/emplar.h"
#include "ff/elec.h"
#include "md.h"
#include "tool/error.h"
#include "tool/zero.h"

namespace tinker {
void emplar(int vers)
{
#if TINKER_CUDART
   bool do_v = vers & calc::virial;

   zeroOnHost(energy_em, virial_em);

   mpole_init(vers);
   emplar_cu(vers);
   torque(vers, demx, demy, demz);
   if (do_v) {
      virial_buffer u2 = vir_trq;
      virial_prec v2[9];
      virial_reduce(v2, u2);
      for (int iv = 0; iv < 9; ++iv)
         virial_elec[iv] += v2[iv];
   }
#else
   (void)vers;
   TINKER_THROW("EMPLAR is only available for PBC systems in CUDA.\n");
#endif
}
}
