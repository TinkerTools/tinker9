#include "ff/hippo/cflux.h"
#include "ff/atom.h"
#include "ff/hippo/elechippo.h"
#include "ff/molecule.h"
#include "ff/nblist.h"
#include "ff/pchg/evalence.h"
#include "ff/potent.h"
#include "tool/io.h"
#include "tool/zero.h"
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
void cfluxData(RcOp op)
{
   if (!potent::use_chgflx)
      return;

   if (op & rc_dealloc) {
      darray::deallocate(bflx, aflx, abflx);
      darray::deallocate(pdelta, atomic, balist);
      darray::deallocate(mono0);
      darray::deallocate(decfx, decfy, decfz, pot);
   }

   if (op & rc_alloc) {
      darray::allocate(n, &pdelta, &atomic);
      darray::allocate(n, &mono0);
      darray::allocate(nbond, &bflx);
      darray::allocate(nangle, &aflx, &abflx, &balist);

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
      darray::copyin(g::q0, nbond, bflx, cflux::bflx);
      darray::copyin(g::q0, nangle, aflx, cflux::aflx);
      darray::copyin(g::q0, nangle, abflx, cflux::abflx);
      darray::copyin(g::q0, n, atomic, atomid::atomic);
      darray::copyin(g::q0, n, mono0, mpole::mono0);

      if (rc_flag & calc::grad)
         darray::zero(g::q0, n, decfx, decfy, decfz, pot);

      std::vector<int> ibalstvec(nangle * 2);
      for (size_t i = 0; i < ibalstvec.size(); ++i) {
         ibalstvec[i] = atmlst::balist[i] - 1;
      }
      darray::copyin(g::q0, nangle, balist, ibalstvec.data());
      waitFor(g::q0);
   }
}

void cfluxZeroPot()
{
   darray::zero(g::q0, n, pot);
}

extern void alterchg_acc();
void alterchg()
{
   alterchg_acc();
}

extern void dcflux_acc(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz, VirialBuffer vir);
void dcflux(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz, VirialBuffer v)
{
   dcflux_acc(vers, gx, gy, gz, v);
}
}
