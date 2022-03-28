#include "add.h"
#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "seq/strtor.h"

namespace tinker {
template <class Ver>
void estrtor_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = bufferSize();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,ebt,vir_ebt,debtx,debty,debtz,\
               ist,kst,bl,itors,tors1,tors2,tors3)
   for (int istrtor = 0; istrtor < nstrtor; ++istrtor) {
      int offset = istrtor & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_strtor<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz, debtx, debty, debtz,

         storunit, istrtor, ist, kst,

         bl, itors, tors1, tors2, tors3,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, ebt, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ebt, offset);
   }
}

void estrtor_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      estrtor_acc1<calc::V0>();
   else if (vers == calc::v1)
      estrtor_acc1<calc::V1>();
   else if (vers == calc::v4)
      estrtor_acc1<calc::V4>();
   else if (vers == calc::v5)
      estrtor_acc1<calc::V5>();
   else if (vers == calc::v6)
      estrtor_acc1<calc::V6>();
}
}
