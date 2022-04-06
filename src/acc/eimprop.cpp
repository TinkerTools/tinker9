#include "ff/energy.h"
#include "ff/evalence.h"
#include "seq/improp.h"

namespace tinker {
template <class Ver>
void eimprop_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = bufferSize();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,deidx,deidy,deidz,\
               iiprop,kprop,vprop,\
               eid,vir_eid)
   for (int i = 0; i < niprop; ++i) {
      int offset = i & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_improp<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         deidx, deidy, deidz,

         idihunit, i, iiprop, kprop, vprop,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, eid, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eid, offset);
   }
}

void eimprop_acc(int vers)
{
   if (vers == calc::v0 or vers == calc::v3)
      eimprop_acc1<calc::V0>();
   else if (vers == calc::v1)
      eimprop_acc1<calc::V1>();
   else if (vers == calc::v4)
      eimprop_acc1<calc::V4>();
   else if (vers == calc::v5)
      eimprop_acc1<calc::V5>();
   else if (vers == calc::v6)
      eimprop_acc1<calc::V6>();
}
}
