#include "add.h"
#include "evalence.h"
#include "md.h"
#include "seq_pitors.h"

namespace tinker {
template <class Ver>
void epitors_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,deptx,depty,deptz,\
               ipit,kpit,\
               ept,vir_ept)
   for (int i = 0; i < npitors; ++i) {
      int offset = i & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_pitors<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         deptx, depty, deptz,

         ptorunit, i, ipit, kpit,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, ept, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ept, offset);
   } // end for (int i)
}

void epitors_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      epitors_acc1<calc::V0>();
   else if (vers == calc::v1)
      epitors_acc1<calc::V1>();
   else if (vers == calc::v4)
      epitors_acc1<calc::V4>();
   else if (vers == calc::v5)
      epitors_acc1<calc::V5>();
   else if (vers == calc::v6)
      epitors_acc1<calc::V6>();
}
}
