#include "accmanaged.h"


namespace tinker {
void accmanaged_data(rc_op op)
{
   using namespace detail;


   if (op & rc_dealloc) {
      #pragma acc exit data async delete(exx,eyy,ezz,exy,eyz,ezx,\
              vtot1,vtot2,vtot3)
   }


   if (op & rc_alloc) {
      #pragma acc enter data async create(exx,eyy,ezz,exy,eyz,ezx,\
              vtot1,vtot2,vtot3)
   }
}


namespace detail {
energy_prec exx, eyy, ezz, exy, eyz, ezx; // kinetic
vel_prec vtot1, vtot2, vtot3;             // mdrest
}
}
