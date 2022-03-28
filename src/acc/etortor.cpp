#include "add.h"
#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "math/inc.h"
#include "seq/tortor.h"

namespace tinker {
template <class Ver>
void etortor_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = bufferSize();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,dettx,detty,dettz,\
               ibitor,itt,chkttor_ia_,\
               tnx,tny,ttx,tty,tbf,tbx,tby,tbxy,\
               ett,vir_ett)
   for (int itortor = 0; itortor < ntortor; ++itortor) {
      int offset = itortor & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_tortor<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         dettx, detty, dettz,

         ttorunit, itortor, itt, ibitor, chkttor_ia_,

         tnx, tny, ttx, tty, tbf, tbx, tby, tbxy,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, ett, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ett, offset);
   } // end for (int itortor)
}

void etortor_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      etortor_acc1<calc::V0>();
   else if (vers == calc::v1)
      etortor_acc1<calc::V1>();
   else if (vers == calc::v4)
      etortor_acc1<calc::V4>();
   else if (vers == calc::v5)
      etortor_acc1<calc::V5>();
   else if (vers == calc::v6)
      etortor_acc1<calc::V6>();
}
}
