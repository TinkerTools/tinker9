#include "cflux.h"
#include "eangle.h"
#include "ebond.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_fort_str.h"
#include <cassert>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/angpot.hh>
#include <tinker/detail/atmlst.hh>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/cflux.hh>
#include <tinker/detail/chgpen.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/mpole.hh>
#include <tinker/detail/potent.hh>


namespace tinker {
void cflux_data(rc_op op)
{
   if (!potent::use_chgflx)
      return;


   if (op & rc_dealloc) {
      darray::deallocate(bflx, aflx, abflx);
      darray::deallocate(pdelta, atomic, balist);
      darray::deallocate(mono0);
      darray::deallocate(couple_i12, couple_n12);
      darray::deallocate(decfx, decfy, decfz, pot);
   }

   if (op & rc_alloc) {
      darray::allocate(n, &pdelta, &atomic);
      darray::allocate(n, &mono0);
      darray::allocate(nbond, &bflx);
      darray::allocate(nangle, &aflx, &abflx, &balist);
      darray::allocate(n, &couple_i12, &couple_n12);

      if (rc_flag & calc::grad) {
         darray::allocate(n, &decfx, &decfy, &decfz, &pot);
      } else {
         decfx = nullptr;
         decfy = nullptr;
         decfz = nullptr;
         pot = nullptr;
      }
   }

   if (op & rc_init) {
      darray::copyin(WAIT_NEW_Q, nbond, bflx, cflux::bflx);
      darray::copyin(WAIT_NEW_Q, nangle, aflx, cflux::aflx);
      darray::copyin(WAIT_NEW_Q, nangle, abflx, cflux::abflx);
      darray::copyin(WAIT_NEW_Q, n, atomic, atomid::atomic);
      darray::copyin(WAIT_NEW_Q, n, mono0, mpole::mono0);


      if (rc_flag & calc::grad)
         darray::zero(PROCEED_NEW_Q, n, decfx, decfy, decfz, pot);

         
      std::vector<int> ibalstvec(nangle * 2);
      for (size_t i = 0; i < ibalstvec.size(); ++i) {
         ibalstvec[i] = atmlst::balist[i] - 1;
      }

      darray::copyin(WAIT_NEW_Q, nangle, balist, ibalstvec.data());
      
      int c_maxn12 = 8;
      std::vector<int> ibuf;
      ibuf.resize(c_maxn12 * n);
      for (int i = 0; i < n; ++i) {
         int nn = couple::n12[i];
         int base = i * c_maxn12;
         for (int j = 0; j < nn; ++j) {
            int k = couple::i12[i][j];
            ibuf[base + j] = k - 1;
         }
         for (int j = nn; j < c_maxn12; ++j) {
            ibuf[base + j] = -1;
         }
      }
      darray::copyin(WAIT_NEW_Q, n, couple_i12, ibuf.data());
      darray::copyin(WAIT_NEW_Q, n, couple_n12, couple::n12);
   }
}


void alterchg()
{
   alterchg_acc();
}

void dcflux(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz,
            virial_buffer v)
{
   dcflux_acc(vers, gx, gy, gz, v);
}

}
