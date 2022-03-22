#include "add.h"
#include "ff/pchg/evalence.h"
#include "md/md.h"
#include "seq/bond.h"
#include <cassert>

namespace tinker {
template <class Ver>
void ebond_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,debx,deby,debz,\
               ibnd,bl,bk,\
               eb,vir_eb)
   for (int i = 0; i < nbond; ++i) {
      int offset = i & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_bond<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         debx, deby, debz,

         bndtyp, bndunit, i, ibnd, bl, bk, cbnd, qbnd,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, eb, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eb, offset);
   } // end for (int i)
}

void ebond_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      ebond_acc1<calc::V0>();
   else if (vers == calc::v1)
      ebond_acc1<calc::V1>();
   else if (vers == calc::v4)
      ebond_acc1<calc::V4>();
   else if (vers == calc::v5)
      ebond_acc1<calc::V5>();
   else if (vers == calc::v6)
      ebond_acc1<calc::V6>();
   else
      assert(false);
}
}
