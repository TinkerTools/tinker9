#include "eangtor.h"
#include "md.h"
#include "potent.h"
#include "tool/host_zero.h"
#include <tinker/detail/angtor.hh>


namespace tinker {
void eangtor_data(rc_op op)
{
   if (not use_potent(angtor_term))
      return;


   bool rc_a = rc_flag & calc::analyz;


   if (op & rc_dealloc) {
   }


   if (op & rc_alloc) {
   }


   if (op & rc_init) {
   }
}


void eangtor(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_eat, virial_eat);
      size_t bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, eat);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eat);
      if (do_g)
         darray::zero(g::q0, n, deatx, deaty, deatz);
   }


   eangtor_acc(vers);


   if (rc_a) {
      if (do_e) {
         energy_eat = energy_reduce(eat);
         energy_valence += energy_eat;
      }
      if (do_v) {
         virial_reduce(virial_eat, vir_eat);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eat[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deatx, deaty, deatz);
   }
}
}
