#include "add.h"
#include "evalence.h"
#include "md.h"
#include "seq_urey.h"

namespace tinker {
template <class Ver>
void eurey_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;

   auto bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,deubx,deuby,deubz,\
               iury,uk,ul,\
               eub,vir_eub)
   for (int i = 0; i < nurey; ++i) {
      int offset = i & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_urey<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         deubx, deuby, deubz,

         ureyunit, i, iury, uk, ul, cury, qury, x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, eub, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eub, offset);
   } // end for (int i)
}

void eurey_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      eurey_acc1<calc::V0>();
   else if (vers == calc::v1)
      eurey_acc1<calc::V1>();
   else if (vers == calc::v4)
      eurey_acc1<calc::V4>();
   else if (vers == calc::v5)
      eurey_acc1<calc::V5>();
   else if (vers == calc::v6)
      eurey_acc1<calc::V6>();
}
}
