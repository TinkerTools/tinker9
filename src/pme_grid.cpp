#include "error.h"
#include "named_struct.h"
#include "platform.h"
#include "pme.h"


TINKER_NAMESPACE_BEGIN
namespace pltfm_acc {
template <class T>
void grid_put(PMEUnit, real*, real*);
}


namespace pltfm_cu {
void grid_mpole(PMEUnit, real (*)[10]);
void grid_uind(PMEUnit, real (*)[3], real (*)[3]);
}


void grid_mpole(PMEUnit pme_u, real (*fmp)[10])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("grid_mpole(): bsorder is {}; must be 5.\n", bso));


   real* opt1 = reinterpret_cast<real*>(fmp);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM) {
      pltfm_cu::grid_mpole(pme_u, fmp);
   } else
#endif
      pltfm_acc::grid_put<MPOLE>(pme_u, opt1, nullptr);
}


void grid_uind(PMEUnit pme_u, real (*fuind)[3], real (*fuinp)[3])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("grid_uind(): bsorder is {}; must be 5.\n", bso));


   real* opt1 = reinterpret_cast<real*>(fuind);
   real* opt2 = reinterpret_cast<real*>(fuinp);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM) {
      pltfm_cu::grid_uind(pme_u, fuind, fuinp);
   } else
#endif
      pltfm_acc::grid_put<UIND>(pme_u, opt1, opt2);
}
TINKER_NAMESPACE_END
