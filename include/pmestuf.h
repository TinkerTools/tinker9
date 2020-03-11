#pragma once
#include "pme.h"


TINKER_NAMESPACE_BEGIN
void bspline_fill(PMEUnit, int level);


void grid_mpole(PMEUnit, real (*fmp)[10]);
void grid_uind(PMEUnit, real (*find)[3], real (*finp)[3]);


void pme_conv(PMEUnit);
void pme_conv(PMEUnit, virial_buffer v);


void fphi_mpole(PMEUnit);
void fphi_uind(PMEUnit, real (*fdip_phi1)[10], real (*fdip_phi2)[10],
               real (*fdip_sum_phi)[20]);
void fphi_uind2(PMEUnit, real (*fdip_phi1)[10], real (*fdip_phi2)[10]);


void rpole_to_cmp();
void cmp_to_fmp(PMEUnit, const real (*cmp)[10], real (*fmp)[10]);
void cuind_to_fuind(PMEUnit, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3]);
void fphi_to_cphi(PMEUnit, const real (*fphi)[20], real (*cphi)[10]);


//====================================================================//


void pme_cuda_func_config();


void bspline_fill_cu(PMEUnit, int level);


#define TINKER_CU_THETA_ON_THE_FLY_GRID_MPOLE 1
#define TINKER_CU_THETA_ON_THE_FLY_GRID_UIND  0


void grid_mpole_acc(PMEUnit, real (*)[10]);
void grid_mpole_cu(PMEUnit, real (*)[10]);
void grid_uind_acc(PMEUnit, real (*)[3], real (*)[3]);
void grid_uind_cu(PMEUnit, real (*)[3], real (*)[3]);


void pme_conv_acc(PMEUnit, virial_buffer);


void fphi_mpole_acc(PMEUnit, real (*)[20]);
void fphi_mpole_cu(PMEUnit, real (*)[20]);
void fphi_uind_acc(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
void fphi_uind_cu(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
void fphi_uind2_acc(PMEUnit, real (*)[10], real (*)[10]);
void fphi_uind2_cu(PMEUnit, real (*)[10], real (*)[10]);


void rpole_to_cmp_acc();
void cmp_to_fmp_acc(PMEUnit, const real (*)[10], real (*)[10]);
void cuind_to_fuind_acc(PMEUnit, const real (*)[3], const real (*)[3],
                        real (*)[3], real (*)[3]);
void fphi_to_cphi_acc(PMEUnit, const real (*)[20], real (*)[10]);
TINKER_NAMESPACE_END
