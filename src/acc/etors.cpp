#include "add.h"
#include "evalence.h"
#include "math/inc.h"
#include "md.h"
#include "seq/torsion.h"

namespace tinker {
template <class Ver>
void etors_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,detx,dety,detz,\
               itors,tors1,tors2,tors3,tors4,tors5,tors6,\
               et,vir_et)
   for (int i = 0; i < ntors; ++i) {
      int offset = i & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_tors<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         detx, dety, detz,

         torsunit, i, itors,

         tors1, tors2, tors3, tors4, tors5, tors6,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, et, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_et, offset);
   } // end for (int i)
}

void etors_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      etors_acc1<calc::V0>();
   else if (vers == calc::v1)
      etors_acc1<calc::V1>();
   else if (vers == calc::v4)
      etors_acc1<calc::V4>();
   else if (vers == calc::v5)
      etors_acc1<calc::V5>();
   else if (vers == calc::v6)
      etors_acc1<calc::V6>();
}
}
