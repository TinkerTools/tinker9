#include "ff/hippo/expol.h"
#include "ff/atom.h"
#include "ff/hippomod.h"
#include "tool/darray.h"
#include <tinker/detail/expol.hh>
#include <tinker/detail/polpot.hh>

namespace tinker {
void expolData(RcOp op)
{
   // TODO Use format like "Potent::EXPOL"
   if (not polpot::use_expol)
      return;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(kpep, prepep, dmppep, lpep);
   }

   if (op & RcOp::ALLOC) {
      darray::allocate(n, &kpep, &prepep, &dmppep, &lpep);
   }

   if (op & RcOp::INIT) {
      darray::copyin(g::q0, n, kpep, expol::kpep);
      darray::copyin(g::q0, n, prepep, expol::prepep);
      darray::copyin(g::q0, n, dmppep, expol::dmppep);
      darray::copyin(g::q0, n, lpep, expol::lpep);
   }
}
}
