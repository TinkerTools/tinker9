#ifndef TINKER_GPU_DECL_PME_H_
#define TINKER_GPU_DECL_PME_H_

#include "mod_pme.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {

namespace detail_ {
std::vector<pme_st>& pme_objs();
std::vector<pme_st*>& pme_deviceptrs();
}

pme_st& pme_obj(int pme_unit);
pme_st* pme_deviceptr(int pme_unit);

/// This function must be called after pme_data has been called because it
/// needs to know the number of pme objects created.
void fft_data(rc_t rc);
void fftfront(int pme_unit);
void fftback(int pme_unit);

int use_ewald();
void pme_init(int vers);
void pme_data(rc_t rc);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
namespace gpu {
/**
 * @brief
 * make the scalar summation over reciprocal lattice
 */
void pme_conv0(int pme_unit);                 // without virial
void pme_conv1(int pme_unit, real* gpu_vir9); // with virial

void rpole_to_cmp();
/**
 * @brief
 * Input: cmp, cartesian rotated mpole.
 * Output: fmp, fractional rotated mpole.
 */
void cmp_to_fmp(int pme_unit, const real (*cmp)[10], real (*fmp)[10]);
void cuind_to_fuind(int pme_unit, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3]);
/**
 * @brief
 * Input: fphi.
 * Output: cphi.
 */
void fphi_to_cphi(int pme_unit, const real (*fphi)[20], real (*cphi)[10]);
/**
 * @brief
 * Input: fmp.
 * Output: qgrid.
 */
void grid_mpole(int pme_unit, real (*gpu_fmp)[10]);
/**
 * @brief
 * Input: qgrid.
 * Output: fphi.
 */
void fphi_mpole(int pme_unit, real (*gpu_fphi)[20]);

void grid_uind(int pme_unit, real (*gpu_find)[3], real (*gpu_finp)[3]);
void fphi_uind(int pme_unit, real (*gpu_fdip_phi1)[10],
               real (*gpu_fdip_phi2)[10], real (*gpu_fdip_sum_phi)[20]);
void fphi_uind2(int pme_unit, real (*gpu_fdip_phi1)[10],
                real (*gpu_fdip_phi2)[10]);
}
TINKER_NAMESPACE_END

#endif
