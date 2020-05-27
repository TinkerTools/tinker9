#include "eopbend.h"
#include "io_fort_str.h"
#include "md.h"
#include "mod.energi.h"
#include "potent.h"
#include <cassert>
#include <tinker/detail/angpot.hh>
#include <tinker/detail/opbend.hh>

namespace tinker {
void eopbend_data(rc_op op)
{
   if (!use_potent(opbend_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(iopb, opbk);

      buffer_deallocate(rc_flag, eopb, vir_eopb);
      buffer_deallocate(rc_flag & ~calc::analyz, deopbx, deopby, deopbz);
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &iopb, &opbk);

      nopbend = count_bonded_term(opbend_term);
      buffer_allocate(rc_flag, &eopb, &vir_eopb);
      buffer_allocate(rc_flag & ~calc::analyz, &deopbx, &deopby, &deopbz);
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
      darray::copyin(WAIT_NEW_Q, nangle, iopb, ibuf.data());
      darray::copyin(WAIT_NEW_Q, nangle, opbk, opbend::opbk);
      opbunit = angpot::opbunit;
      copb = angpot::copb;
      qopb = angpot::qopb;
      popb = angpot::popb;
      sopb = angpot::sopb;
   }
}

void eopbend(int vers)
{
   eopbend_acc(vers);


   if (rc_flag & calc::analyz) {
      if (vers & calc::energy) {
         energy_eopb = energy_reduce(eopb);
         energy_valence += energy_eopb;
      }
      if (vers & calc::virial) {
         virial_reduce(virial_eopb, vir_eopb);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eopb[iv];
      }
   }
   if (vers & calc::analyz)
      if (vers & calc::grad)
         sum_gradient(gx, gy, gz, deopbx, deopby, deopbz);
}
}
