#ifndef TINKER_GPU_DECL_PME_H_
#define TINKER_GPU_DECL_PME_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/**
 * pme.f and kewald.f
 *
 * allocate igrid(3,n)                   !! integer
 * allocate (bsmod1(nfft1))              !! real*8
 * allocate (bsmod2(nfft2))              !! real*8
 * allocate (bsmod3(nfft3))              !! real*8
 * allocate (bsbuild(bsorder,bsorder))   !! real*8
 * allocate (thetai1(4,bsorder,n))       !! real*8
 * allocate (thetai2(4,bsorder,n))       !! real*8
 * allocate (thetai3(4,bsorder,n))       !! real*8
 * allocate (qgrid(2,nfft1,nfft2,nfft3)) !! real*8
 * allocate (qfac(nfft1,nfft2,nfft3))    !! real*8
 */
struct pme_st {
  real aewald;
  int nfft1, nfft2, nfft3, bsorder;
  int* igrid;                     // deviceptr
  real *bsmod1, *bsmod2, *bsmod3; // deviceptr
  real* qgrid;                    // deviceptr
};

namespace detail_ {
std::vector<pme_st>& pme_objs();
std::vector<pme_st*>& pme_deviceptrs();
}

int pme_open_unit(double aewald, int nfft1, int nfft2, int nfft3, int bsorder);
pme_st& pme_obj(int pme_unit);
pme_st* pme_deviceptr(int pme_unit);

namespace detail_ {
/// This function must be called after pme_data has been called because it
/// needs to know the number of pme objects created.
void fft_data(int op);
}
void fftfront(int pme_unit);
void fftback(int pme_unit);

extern int epme_unit; // electrostatic
extern int ppme_unit; // polarization
extern int dpme_unit; // dispersion

extern double ewald_switch_cut, ewald_switch_off;

extern real (*fmp)[10];
extern real (*cphi)[10];
extern real (*fphi)[20];
void pme_data(int op);
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

/**
 * @brief
 * Input: cmp, cartesian rotated mpole.
 * Output: fmp, fractional rotated mpole.
 */
void cmp_to_fmp(int pme_unit, real (*fmp)[10]);
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
}
TINKER_NAMESPACE_END

#endif
