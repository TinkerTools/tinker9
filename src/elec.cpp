#include "ff/elec.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/energy.h"
#include "ff/pme.h"

namespace tinker {
void chkpole_acc();
void rotpole_acc();
void torque_acc(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz);
static void chkpole()
{
   chkpole_acc();
}

static void rotpole()
{
   rotpole_acc();
}

void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz)
{
   // #if TINKER_CUDART
   //    if (pltfm_config & Platform::CUDA)
   //    else
   // #endif
   torque_acc(vers, dx, dy, dz);
}

void mpoleInit(int vers)
{
   if (vers & calc::grad)
      darray::zero(g::q0, n, trqx, trqy, trqz);
   if (vers & calc::virial)
      darray::zero(g::q0, bufferSize(), vir_trq);

   chkpole();
   rotpole();

   if (useEwald()) {
      rpole_to_cmp();
      if (vir_m)
         darray::zero(g::q0, bufferSize(), vir_m);
      if (pltfm_config & Platform::CUDA) {
         bool precompute_theta =
            (!TINKER_CU_THETA_ON_THE_FLY_GRID_MPOLE) || (!TINKER_CU_THETA_ON_THE_FLY_GRID_UIND);
         if (epme_unit.valid()) {
            if (precompute_theta)
               bspline_fill(epme_unit, 3);
         }
         if (ppme_unit.valid() && (ppme_unit != epme_unit)) {
            if (precompute_theta)
               bspline_fill(ppme_unit, 2);
         }
         if (pvpme_unit.valid()) {
            if (precompute_theta)
               bspline_fill(pvpme_unit, 2);
         }
      }
   }
}
}
