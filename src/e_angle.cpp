#include "e_angle.h"
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
      device_array::deallocate(iang, ak, anat, angtyp);

      buffer_deallocate(ea, vir_ea);
   }

   if (op & rc_alloc) {
      nangle = count_bonded_term(angle_term);
      device_array::allocate(nangle, &iang, &ak, &anat, &angtyp);

      buffer_allocate(&ea, &vir_ea);
   }

   if (op & rc_init) {
      std::vector<int> iangvec(nangle * 4);
      for (size_t i = 0; i < iangvec.size(); ++i) {
         iangvec[i] = angbnd::iang[i] - 1;
      }
      device_array::copyin(nangle, iang, iangvec.data());
      device_array::copyin(nangle, ak, angbnd::ak);
      device_array::copyin(nangle, anat, angbnd::anat);

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
      device_array::copyin(nangle, angtyp, angtypvec.data());
   }
}

void eangle(int vers)
{
   extern void eangle_acc(int);
   eangle_acc(vers);
}
TINKER_NAMESPACE_END
