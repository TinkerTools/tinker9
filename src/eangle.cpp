#include "eangle.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include <cassert>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/angpot.hh>

TINKER_NAMESPACE_BEGIN
void eangle_data(rc_op op)
{
   if (!use_potent(angle_term) && !use_potent(strbnd_term) &&
       !use_potent(opbend_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(iang, ak, anat, angtyp);

      buffer_deallocate(rc_flag, ea, deax, deay, deaz, vir_ea);
   }

   if (op & rc_alloc) {
      nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &iang, &ak, &anat, &angtyp);

      buffer_allocate(rc_flag, &ea, &deax, &deay, &deaz, &vir_ea, &energy_ea);
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
}
TINKER_NAMESPACE_END
