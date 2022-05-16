#include "ff/amoeba/empole.h"
#include "ff/amoebamod.h"
#include "ff/energy.h"
#include "math/zero.h"
#include "tool/externfunc.h"

namespace tinker {
TINKER_FVOID2(acc0, cu1, emplar, int);
void emplar(int vers)
{
   bool do_v = vers & calc::virial;

   zeroOnHost(energy_em, virial_em);

   mpoleInit(vers);
   TINKER_FCALL2(acc0, cu1, emplar, vers);
   torque(vers, demx, demy, demz);
   if (do_v) {
      VirialBuffer u2 = vir_trq;
      virial_prec v2[9];
      virialReduce(v2, u2);
      for (int iv = 0; iv < 9; ++iv)
         virial_elec[iv] += v2[iv];
   }
}
}
