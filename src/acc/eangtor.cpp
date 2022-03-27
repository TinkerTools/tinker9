#include "add.h"
#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "seq/angtor.h"

namespace tinker {
template <class Ver>
void eangtor_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,eat,vir_eat,deatx,deaty,deatz,\
               iat,kant,anat,itors,tors1,tors2,tors3)
   for (int iangtor = 0; iangtor < nangtor; ++iangtor) {
      int offset = iangtor & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_angtor<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz, deatx, deaty, deatz,

         atorunit, iangtor, iat, kant,

         anat, itors, tors1, tors2, tors3,

         x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, eat, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eat, offset);
   }
}

void eangtor_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      eangtor_acc1<calc::V0>();
   else if (vers == calc::v1)
      eangtor_acc1<calc::V1>();
   else if (vers == calc::v4)
      eangtor_acc1<calc::V4>();
   else if (vers == calc::v5)
      eangtor_acc1<calc::V5>();
   else if (vers == calc::v6)
      eangtor_acc1<calc::V6>();
}
}
