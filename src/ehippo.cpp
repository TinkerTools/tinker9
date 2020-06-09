#include "ehippo.h"
#include "nblist.h"


namespace tinker {
void ehippo(int vers)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      ehippo_cu(vers);


   if (vers & calc::energy) {
      energy_buffer u = ect;
      energy_ect = energy_reduce(u);
      energy_elec += energy_ect;
   }
   if (vers & calc::virial) {
      virial_buffer u = vir_ect;
      virial_prec v1[9];
      virial_reduce(v1, u);
      for (int iv = 0; iv < 9; ++iv) {
         virial_elec[iv] += v1[iv];
      }
   }
#endif
}
}
