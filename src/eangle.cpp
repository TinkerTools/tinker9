#include "eangle.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include <cassert>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/angpot.hh>

namespace tinker {
void eangle_data(rc_op op)
{
   if (!use_potent(angle_term) && !use_potent(strbnd_term) &&
       !use_potent(opbend_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(iang, ak, anat, angtyp);

      buffer_deallocate(rc_flag, ea, vir_ea);
      buffer_deallocate(rc_flag & ~calc::analyz, deax, deay, deaz);
   }

   if (op & rc_alloc) {
      nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &iang, &ak, &anat, &angtyp);

      buffer_allocate(rc_flag, &ea, &vir_ea);
      buffer_allocate(rc_flag & ~calc::analyz, &deax, &deay, &deaz);
   }

   if (op & rc_init) {
      std::vector<int> iangvec(nangle * 4);
      for (size_t i = 0; i < iangvec.size(); ++i) {
         iangvec[i] = angbnd::iang[i] - 1;
      }
      darray::copyin(WAIT_NEW_Q, nangle, iang, iangvec.data());
      darray::copyin(WAIT_NEW_Q, nangle, ak, angbnd::ak);
      darray::copyin(WAIT_NEW_Q, nangle, anat, angbnd::anat);

      angunit = angpot::angunit;
      cang = angpot::cang;
      qang = angpot::qang;
      pang = angpot::pang;
      sang = angpot::sang;
      std::vector<eangle_t> angtypvec(nangle);
      for (int i = 0; i < nangle; ++i) {
         fstr_view atyp = angpot::angtyp[i];
         if (atyp == "IN-PLANE")
            angtypvec[i] = eangle_t::in_plane;
         else if (atyp == "HARMONIC")
            angtypvec[i] = eangle_t::harmonic;
         else if (atyp == "LINEAR")
            angtypvec[i] = eangle_t::linear;
         else if (atyp == "FOURIER")
            angtypvec[i] = eangle_t::fourier;
         else {
            assert(false);
         }
      }
      darray::copyin(WAIT_NEW_Q, nangle, angtyp, angtypvec.data());
   }
}

void eangle(int vers)
{
   eangle_acc(vers);


   if (rc_flag & calc::analyz) {
      if (vers & calc::energy) {
         energy_ea = energy_reduce(ea);
         energy_valence += energy_ea;
      }
      if (vers & calc::virial) {
         virial_reduce(virial_ea, vir_ea);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ea[iv];
      }
   }
   if (vers & calc::analyz)
      if (vers & calc::grad)
         sum_gradient(gx, gy, gz, deax, deay, deaz);
}
}
