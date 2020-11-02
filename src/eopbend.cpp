#include "eopbend.h"
#include "md.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_fort_str.h"
#include <cassert>
#include <tinker/detail/angpot.hh>
#include <tinker/detail/opbend.hh>

namespace tinker {
void eopbend_data(rc_op op)
{
   if (!use_potent(opbend_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(iopb, opbk);

      if (rc_a)
         buffer_deallocate(rc_flag, eopb, vir_eopb, deopbx, deopby, deopbz);
      eopb = nullptr;
      vir_eopb = nullptr;
      deopbx = nullptr;
      deopby = nullptr;
      deopbz = nullptr;
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &iopb, &opbk);
      nopbend = count_bonded_term(opbend_term);

      eopb = eng_buf;
      vir_eopb = vir_buf;
      deopbx = gx;
      deopby = gy;
      deopbz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &eopb, &vir_eopb, &deopbx, &deopby, &deopbz);
   }

   if (op & rc_init) {
      fstr_view otyp = angpot::opbtyp;
      if (otyp == "W-D-C")
         opbtyp = eopbend_t::w_d_c;
      else if (otyp == "ALLINGER")
         opbtyp = eopbend_t::allinger;
      else
         assert(false);
      int nangle = count_bonded_term(angle_term);
      std::vector<int> ibuf(nangle);
      for (int i = 0; i < nangle; ++i)
         ibuf[i] = opbend::iopb[i] - 1;
      darray::copyin(asyncq, nangle, iopb, ibuf.data());
      darray::copyin(asyncq, nangle, opbk, opbend::opbk);
      wait_for(asyncq);
      opbunit = angpot::opbunit;
      copb = angpot::copb;
      qopb = angpot::qopb;
      popb = angpot::popb;
      sopb = angpot::sopb;
   }
}

void eopbend(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_eopb, virial_eopb);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(asyncq, bsize, eopb);
      if (do_v)
         darray::zero(asyncq, bsize, vir_eopb);
      if (do_g)
         darray::zero(asyncq, n, deopbx, deopby, deopbz);
   }


   eopbend_acc(vers);


   if (rc_a) {
      if (do_e) {
         energy_eopb = energy_reduce(eopb);
         energy_valence += energy_eopb;
      }
      if (do_v) {
         virial_reduce(virial_eopb, vir_eopb);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eopb[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deopbx, deopby, deopbz);
   }
}
}
