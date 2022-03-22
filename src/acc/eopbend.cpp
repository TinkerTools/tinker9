#include "add.h"
#include "evalence.h"
#include "math/inc.h"
#include "md.h"
#include "seq/opbend.h"
#include <cassert>

namespace tinker {
template <class Ver>
void eopbend_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,deopbx,deopby,deopbz,\
               iopb,opbk,iang,\
               eopb,vir_eopb)
   for (int iopbend = 0; iopbend < nopbend; ++iopbend) {
      int offset = iopbend & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_opbend<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         deopbx, deopby, deopbz,

         opbtyp, opbunit, iopbend, iopb, opbk, iang, copb, qopb, popb, sopb,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, eopb, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eopb, offset);
   } // end for (int iopbend)
}

void eopbend_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      eopbend_acc1<calc::V0>();
   else if (vers == calc::v1)
      eopbend_acc1<calc::V1>();
   else if (vers == calc::v4)
      eopbend_acc1<calc::V4>();
   else if (vers == calc::v5)
      eopbend_acc1<calc::V5>();
   else if (vers == calc::v6)
      eopbend_acc1<calc::V6>();
}
}
