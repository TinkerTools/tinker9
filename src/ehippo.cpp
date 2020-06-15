#include "ehippo.h"
#include "md.h"
#include "nblist.h"
#include "tool/host_zero.h"

namespace tinker {
void ehippo(int vers)
{
#if TINKER_CUDART
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   host_zero(energy_ect, virial_ect);
   size_t bsize = buffer_size();
   if (do_a)
      darray::zero(PROCEED_NEW_Q, bsize, nct);
   if (do_e)
      darray::zero(PROCEED_NEW_Q, bsize, ect);
   if (do_v)
      darray::zero(PROCEED_NEW_Q, bsize, vir_ect);
   if (do_g)
      darray::zero(PROCEED_NEW_Q, n, dectx, decty, dectz);


   if (mlist_version() & NBL_SPATIAL)
      ehippo_cu(vers);


   if (do_e) {
      energy_buffer u = ect;
      energy_ect = energy_reduce(u);
      energy_elec += energy_ect;
   }
   if (do_v) {
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
