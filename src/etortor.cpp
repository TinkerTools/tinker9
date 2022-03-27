#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "ff/potent.h"
#include "tool/zero.h"
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bitor.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/torpot.hh>
#include <tinker/detail/tortor.hh>

namespace tinker {
void etortorData(RcOp op)
{
   if (!use_potent(tortor_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(ibitor, itt, tnx, tny, ttx, tty, tbf, tbx, tby, tbxy, chkttor_ia_);

      if (rc_a)
         buffer_deallocate(rc_flag, ett, vir_ett, dettx, detty, dettz);
      ett = nullptr;
      vir_ett = nullptr;
      dettx = nullptr;
      detty = nullptr;
      dettz = nullptr;
   }

   if (op & rc_alloc) {
      nbitor = bitor_::nbitor;
      darray::allocate(nbitor, &ibitor, &itt);

      darray::allocate(ktrtor::maxntt, &tnx, &tny, &ttx, &tty, &tbf, &tbx, &tby, &tbxy);

      ntortor = count_bonded_term(tortor_term);
      darray::allocate(ntortor, &chkttor_ia_);

      ett = eng_buf;
      vir_ett = vir_buf;
      dettx = gx;
      detty = gy;
      dettz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &ett, &vir_ett, &dettx, &detty, &dettz);
   }

   if (op & rc_init) {
      std::vector<int> ibuf;

      ibuf.resize(5 * nbitor);
      for (int i = 0; i < 5 * nbitor; ++i)
         ibuf[i] = bitor_::ibitor[i] - 1;
      darray::copyin(g::q0, nbitor, ibitor, ibuf.data());
      wait_for(g::q0);

      ibuf.resize(3 * nbitor);
      for (int i = 0; i < 3 * nbitor; ++i)
         ibuf[i] = tortor::itt[i] - 1;
      darray::copyin(g::q0, nbitor, itt, ibuf.data());
      wait_for(g::q0);

      ibuf.resize(ktrtor::maxntt);
      for (int i = 0; i < ktrtor::maxntt; ++i)
         ibuf[i] = ktrtor::tnx[i];
      darray::copyin(g::q0, ktrtor::maxntt, tnx, ibuf.data());
      wait_for(g::q0);
      for (int i = 0; i < ktrtor::maxntt; ++i)
         ibuf[i] = ktrtor::tny[i];
      darray::copyin(g::q0, ktrtor::maxntt, tny, ibuf.data());
      darray::copyin(g::q0, ktrtor::maxntt, ttx, &ktrtor::ttx[0][0]);
      darray::copyin(g::q0, ktrtor::maxntt, tty, &ktrtor::tty[0][0]);
      darray::copyin(g::q0, ktrtor::maxntt, tbf, &ktrtor::tbf[0][0]);
      darray::copyin(g::q0, ktrtor::maxntt, tbx, &ktrtor::tbx[0][0]);
      darray::copyin(g::q0, ktrtor::maxntt, tby, &ktrtor::tby[0][0]);
      darray::copyin(g::q0, ktrtor::maxntt, tbxy, &ktrtor::tbxy[0][0]);
      wait_for(g::q0);

      ttorunit = torpot::ttorunit;

      // see also subroutine chkttor in etortor.f
      ibuf.resize(ntortor);
      for (int itortor = 0; itortor < ntortor; ++itortor) {
         int ib, ic, id;
         int i = tortor::itt[itortor * 3] - 1;
         int flag = tortor::itt[itortor * 3 + 2];
         if (flag == 1) {
            ib = bitor_::ibitor[5 * i + 1] - 1;
            ic = bitor_::ibitor[5 * i + 2] - 1;
            id = bitor_::ibitor[5 * i + 3] - 1;
         } else {
            ib = bitor_::ibitor[5 * i + 3] - 1;
            ic = bitor_::ibitor[5 * i + 2] - 1;
            id = bitor_::ibitor[5 * i + 1] - 1;
         }

         int ia = -1; // -1 is the default ia value for non-chiral center ic
         if (couple::n12[ic] == 4) {
            /*
             *  j    k
             *   \  /
             *    ic
             *   /  \
             * ib    id
             */
            int j = -1;
            int k;
            for (int i = 0; i < 4; ++i) {
               int m = couple::i12[ic][i] - 1;
               if (m != ib && m != id) {
                  if (j == -1) // if j is uninitialized
                     j = m;
                  else // j has been uninitialized
                     k = m;
               }
            }

            if (atoms::type[j] > atoms::type[k])
               ia = j;
            if (atoms::type[k] > atoms::type[j])
               ia = k;
            // else ia is uninitialized
            if (atomid::atomic[j] > atomid::atomic[k])
               ia = j;
            if (atomid::atomic[k] > atomid::atomic[j])
               ia = k;
            // else ia is unchanged
         }
         ibuf[itortor] = ia;
      }
      darray::copyin(g::q0, ntortor, chkttor_ia_, ibuf.data());
      wait_for(g::q0);
   }
}

void etortor(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_ett, virial_ett);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, ett);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ett);
      if (do_g)
         darray::zero(g::q0, n, dettx, detty, dettz);
   }

   etortor_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_ett = energy_reduce(ett);
         energy_valence += energy_ett;
      }
      if (do_v) {
         virial_reduce(virial_ett, vir_ett);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ett[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, dettx, detty, dettz);
   }
}
}
