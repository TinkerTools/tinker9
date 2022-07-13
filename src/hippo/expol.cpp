#include "ff/hippo/expol.h"
#include "ff/atom.h"
#include "ff/hippomod.h"
#include "tool/darray.h"
#include <tinker/detail/kexpl.hh>
#include <tinker/detail/polpot.hh>

namespace tinker {
void expolData(RcOp op)
{
   // TODO Return if expol not used
   if (not polpot::use_expol)
      return;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(pepk, peppre, pepdmp, pepl);
   }

   if (op & RcOp::ALLOC) {
      darray::allocate(n, &pepk, &peppre, &pepdmp, &pepl);
   }

   if (op & RcOp::INIT) {
      darray::copyin(g::q0, n, pepk, kexpl::pepk);
      darray::copyin(g::q0, n, peppre, kexpl::peppre);
      darray::copyin(g::q0, n, pepdmp, kexpl::pepdmp);
      darray::copyin(g::q0, n, pepl, kexpl::pepl);
   }
}
}
