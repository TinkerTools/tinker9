#include "ff/hippo/expol.h"
#include "ff/atom.h"
#include "ff/hippomod.h"
#include "tool/darray.h"
#include "tool/iofortstr.h"
#include <tinker/detail/expol.hh>
#include <tinker/detail/polpot.hh>

namespace tinker {
void expolData(RcOp op)
{
   if (not polpot::use_expol)
      return;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(kpep, prepep, dmppep, lpep);
      darray::deallocate(polscale, polinv);

      scrtyp = ExpolScr::NONE;
   }

   if (op & RcOp::ALLOC) {
      darray::allocate(n, &kpep, &prepep, &dmppep, &lpep);
      darray::allocate(n, &polscale, &polinv);
   }

   if (op & RcOp::INIT) {
      darray::copyin(g::q0, n, kpep, expol::kpep);
      darray::copyin(g::q0, n, prepep, expol::prepep);
      darray::copyin(g::q0, n, dmppep, expol::dmppep);
      darray::copyin(g::q0, n, lpep, expol::lpep);

      FstrView scrview = polpot::scrtyp;
      if (scrview == "S2U")
         scrtyp = ExpolScr::S2U;
      else if (scrview == "S2 ")
         scrtyp = ExpolScr::S2;
      else if (scrview == "G  ")
         scrtyp = ExpolScr::G;
      else
         scrtyp = ExpolScr::NONE;
   }
}
}
