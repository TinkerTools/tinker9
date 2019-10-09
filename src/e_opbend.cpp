#include "e_opbend.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include <cassert>
#include <ext/tinker/detail/angpot.hh>
#include <ext/tinker/detail/opbend.hh>

TINKER_NAMESPACE_BEGIN
void eopbend_data(rc_op op)
{
   if (!use_potent(opbend_term))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate(iopb, opbk);

      buffer_deallocate(eopb, vir_eopb);
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      device_array::allocate(nangle, &iopb, &opbk);

      nopbend = count_bonded_term(opbend_term);
      buffer_allocate(&eopb, &vir_eopb);
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
      device_array::copyin(nangle, iopb, ibuf.data());
      device_array::copyin(nangle, opbk, opbend::opbk);
      opbunit = angpot::opbunit;
      copb = angpot::copb;
      qopb = angpot::qopb;
      popb = angpot::popb;
      sopb = angpot::sopb;
   }
}

extern void eopbend_acc_impl_(int vers);
void eopbend(int vers)
{
   eopbend_acc_impl_(vers);
}
TINKER_NAMESPACE_END
