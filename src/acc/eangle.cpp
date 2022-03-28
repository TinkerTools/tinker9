#include "add.h"
#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "math/inc.h"
#include "seq/angle.h"

namespace tinker {
template <class Ver>
void eangle_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = bufferSize();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,deax,deay,deaz,\
               iang,anat,ak,afld,angtyp,\
               ea,vir_ea)
   for (int i = 0; i < nangle; ++i) {
      int offset = i & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_angle<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         deax, deay, deaz,

         angtyp, angunit, i, iang, anat, ak, afld,

         cang, qang, pang, sang,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, ea, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ea, offset);
   } // end for (int i)
}

void eangle_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      eangle_acc1<calc::V0>();
   else if (vers == calc::v1)
      eangle_acc1<calc::V1>();
   else if (vers == calc::v4)
      eangle_acc1<calc::V4>();
   else if (vers == calc::v5)
      eangle_acc1<calc::V5>();
   else if (vers == calc::v6)
      eangle_acc1<calc::V6>();
}
}
