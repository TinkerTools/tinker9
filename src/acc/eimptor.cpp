#include "eimptor.h"
#include "add.h"
#include "md.h"
#include "seq_imptor.h"

namespace tinker {
template <class Ver>
void eimptor_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   size_t bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,deitx,deity,deitz,\
               iitors,itors1,itors2,itors3,\
               eit,vir_eit)
   for (int i = 0; i < nitors; ++i) {
      int offset = i & (bufsize - 1);
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_imptor<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

         deitx, deity, deitz,

         itorunit, i, iitors, itors1, itors2, itors3, x, y, z);
      if CONSTEXPR (do_e)
         atomic_add(e, eit, offset);
      if CONSTEXPR (do_v)
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eit, offset);
   }
}

void eimptor_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      eimptor_acc1<calc::V0>();
   else if (vers == calc::v1)
      eimptor_acc1<calc::V1>();
   else if (vers == calc::v4)
      eimptor_acc1<calc::V4>();
   else if (vers == calc::v5)
      eimptor_acc1<calc::V5>();
   else if (vers == calc::v6)
      eimptor_acc1<calc::V6>();
}
}
