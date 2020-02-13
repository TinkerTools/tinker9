#include "e_tortor.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bitor.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/torpot.hh>
#include <tinker/detail/tortor.hh>

TINKER_NAMESPACE_BEGIN
void etortor_data(rc_op op)
{
   if (!use_potent(tortor_term))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate(ibitor, itt, tnx, tny, ttx, tty, tbf, tbx, tby,
                               tbxy, chkttor_ia_);

      buffer_deallocate(ett, vir_ett);
   }

   if (op & rc_alloc) {
      nbitor = bitor_::nbitor;
      device_array::allocate(nbitor, &ibitor, &itt);

      device_array::allocate(ktrtor::maxntt, &tnx, &tny, &ttx, &tty, &tbf, &tbx,
                             &tby, &tbxy);

      ntortor = count_bonded_term(tortor_term);
      device_array::allocate(ntortor, &chkttor_ia_);

      buffer_allocate(&ett, &vir_ett);
   }

   if (op & rc_init) {
      std::vector<int> ibuf;

      ibuf.resize(5 * nbitor);
      for (int i = 0; i < 5 * nbitor; ++i)
         ibuf[i] = bitor_::ibitor[i] - 1;
      device_array::copyin(WAIT_NEW_Q, nbitor, ibitor, ibuf.data());

      ibuf.resize(3 * nbitor);
      for (int i = 0; i < 3 * nbitor; ++i)
         ibuf[i] = tortor::itt[i] - 1;
      device_array::copyin(WAIT_NEW_Q, nbitor, itt, ibuf.data());

      ibuf.resize(ktrtor::maxntt);
      for (int i = 0; i < ktrtor::maxntt; ++i)
         ibuf[i] = ktrtor::tnx[i];
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, tnx, ibuf.data());
      for (int i = 0; i < ktrtor::maxntt; ++i)
         ibuf[i] = ktrtor::tny[i];
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, tny, ibuf.data());
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, ttx, &ktrtor::ttx[0][0]);
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, tty, &ktrtor::tty[0][0]);
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, tbf, &ktrtor::tbf[0][0]);
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, tbx, &ktrtor::tbx[0][0]);
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, tby, &ktrtor::tby[0][0]);
      device_array::copyin(WAIT_NEW_Q, ktrtor::maxntt, tbxy,
                           &ktrtor::tbxy[0][0]);

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
      device_array::copyin(WAIT_NEW_Q, ntortor, chkttor_ia_, ibuf.data());
   }
}

void etortor(int vers)
{
   extern void etortor_acc(int);
   etortor_acc(vers);
}
TINKER_NAMESPACE_END
