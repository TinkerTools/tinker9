#include "pmestuf.h"
#include "tool/error.h"

namespace tinker {
void bspline_fill(PMEUnit pme_u, int level)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      bspline_fill_cu(pme_u, level);
#else
   (void)pme_u;
   (void)level;
#endif
}


//====================================================================//


void grid_pchg(PMEUnit pme_u, real* pchg)
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("grid_pchg(): bsorder is %d; must be 5.\n", bso));


#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      grid_pchg_cu(pme_u, pchg);
   else
#endif
      grid_pchg_acc(pme_u, pchg);
}


void grid_mpole(PMEUnit pme_u, real (*fmp)[10])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("grid_mpole(): bsorder is %d; must be 5.\n", bso));


#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      grid_mpole_cu(pme_u, fmp);
   else
#endif
      grid_mpole_acc(pme_u, fmp);
}


void grid_uind(PMEUnit pme_u, real (*fuind)[3], real (*fuinp)[3])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("grid_uind(): bsorder is %d; must be 5.\n", bso));


#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      grid_uind_cu(pme_u, fuind, fuinp);
   else
#endif
      grid_uind_acc(pme_u, fuind, fuinp);
}


//====================================================================//


void pme_conv(PMEUnit pme_u)
{
   pme_conv_acc(pme_u, nullptr, nullptr);
}


void pme_conv(PMEUnit pme_u, virial_buffer gpu_vir)
{
   pme_conv_acc(pme_u, nullptr, gpu_vir);
}


void pme_conv(PMEUnit pme_u, energy_buffer gpu_e)
{
   pme_conv_acc(pme_u, gpu_e, nullptr);
}


void pme_conv(PMEUnit pme_u, energy_buffer gpu_e, virial_buffer gpu_vir)
{
   pme_conv_acc(pme_u, gpu_e, gpu_vir);
}


//====================================================================//


void fphi_mpole(PMEUnit pme_u)
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("fphi_mpole(): bsorder is %d; must be 5.\n", bso));


#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      fphi_mpole_cu(pme_u, fphi);
   else
#endif
      fphi_mpole_acc(pme_u, fphi);
}


void fphi_uind(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10],
               real (*fdip_sum_phi)[20])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("fphi_uind(): bsorder is %d; must be 5.\n", bso));


#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      fphi_uind_cu(pme_u, fdip_phi1, fdip_phi2, fdip_sum_phi);
   else
#endif
      fphi_uind_acc(pme_u, fdip_phi1, fdip_phi2, fdip_sum_phi);
}


void fphi_uind2(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10])
{
   int bso = pme_u->bsorder;
   if (bso != 5)
      TINKER_THROW(format("fphi_uind2(): bsorder is %d; must be 5.\n", bso));


#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      fphi_uind2_cu(pme_u, fdip_phi1, fdip_phi2);
   else
#endif
      fphi_uind2_acc(pme_u, fdip_phi1, fdip_phi2);
}


//====================================================================//


void rpole_to_cmp()
{
   rpole_to_cmp_acc();
}


void cmp_to_fmp(PMEUnit pme_u, const real (*cmp)[10], real (*fmp)[10])
{
   cmp_to_fmp_acc(pme_u, cmp, fmp);
}


void cuind_to_fuind(PMEUnit pme_u, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3])
{
   cuind_to_fuind_acc(pme_u, cind, cinp, fuind, fuinp);
}


void fphi_to_cphi(PMEUnit pme_u, const real (*fphi)[20], real (*cphi)[10])
{
   fphi_to_cphi_acc(pme_u, fphi, cphi);
}
}
