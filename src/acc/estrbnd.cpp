#include "ff/energy.h"
#include "ff/evalence.h"
#include "seq/strbnd.h"

namespace tinker {
template <class Ver>
void estrbnd_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = bufferSize();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,debax,debay,debaz,\
               isb,sbk,iang,anat,bl,\
               eba,vir_eba)
   for (int istrbnd = 0; istrbnd < nstrbnd; ++istrbnd) {
      int offset = istrbnd & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_strbnd<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         debax, debay, debaz,

         stbnunit, istrbnd, isb, sbk,

         bl, iang, anat,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, eba, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eba, offset);
   }
}

void estrbnd_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      estrbnd_acc1<calc::V0>();
   else if (vers == calc::v1)
      estrbnd_acc1<calc::V1>();
   else if (vers == calc::v4)
      estrbnd_acc1<calc::V4>();
   else if (vers == calc::v5)
      estrbnd_acc1<calc::V5>();
   else if (vers == calc::v6)
      estrbnd_acc1<calc::V6>();
}
}
