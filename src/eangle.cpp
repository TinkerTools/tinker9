#include "eangle.h"
#include "md.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_fort_str.h"
#include <cassert>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/angpot.hh>
#include <tinker/detail/potent.hh>

namespace tinker {
void eangle_data(rc_op op)
{
   if (not use_potent(angle_term) and not use_potent(strbnd_term) and
       not use_potent(opbend_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(iang, ak, anat, afld, angtyp);

      if (rc_a)
         buffer_deallocate(rc_flag, ea, vir_ea, deax, deay, deaz);
      ea = nullptr;
      vir_ea = nullptr;
      deax = nullptr;
      deay = nullptr;
      deaz = nullptr;
   }

   if (op & rc_alloc) {
      nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &iang, &ak, &anat, &afld, &angtyp);

      ea = eng_buf;
      vir_ea = vir_buf;
      deax = gx;
      deay = gy;
      deaz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &ea, &vir_ea, &deax, &deay, &deaz);
   }

   if (op & rc_init) {
      std::vector<int> iangvec(nangle * 4);
      for (size_t i = 0; i < iangvec.size(); ++i) {
         iangvec[i] = angbnd::iang[i] - 1;
      }
      darray::copyin(g::q0, nangle, iang, iangvec.data());
      darray::copyin(g::q0, nangle, ak, angbnd::ak);
      darray::copyin(g::q0, nangle, anat, angbnd::anat);
      darray::copyin(g::q0, nangle, afld, angbnd::afld);
      wait_for(g::q0);

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
      darray::copyin(g::q0, nangle, angtyp, angtypvec.data());
      wait_for(g::q0);
   }
}

void eangle(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_ea, virial_ea);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, ea);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ea);
      if (do_g)
         darray::zero(g::q0, n, deax, deay, deaz);
   }


   eangle_acc(vers);


   if (rc_a) {
      if (do_e) {
         energy_ea = energy_reduce(ea);
         energy_valence += energy_ea;
      }
      if (do_v) {
         virial_reduce(virial_ea, vir_ea);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ea[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deax, deay, deaz);
   }
}
}
