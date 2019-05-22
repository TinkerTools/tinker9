#ifndef TINKER_GPU_DECL_PME_H_
#define TINKER_GPU_DECL_PME_H_

#include "decl_real.h"
#include <vector>

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
  int nfft1, nfft2, nfft3, bsorder;
  int* igrid;                                       // deviceptr
  real *bsmod1, *bsmod2, *bsmod3, *bsbuild;         // deviceptr
  real *thetai1, *thetai2, *thetai3, *qgrid, *qfac; // deviceptr
};

namespace detail_ {
std::vector<pme_st>& pme_objs();
std::vector<pme_st*>& pme_deviceptrs();
}

int pme_open_unit(int nfft1, int nfft2, int nfft3, int bsorder);
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

void pme_data(int op);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
namespace gpu {
void cmp_to_fmp(const real (*cmp)[10], real (*fmp)[10], int pme_unit);
void fphi_to_cphi(const real (*fphi)[20], real (*cphi)[10], int pme_unit);
}
TINKER_NAMESPACE_END

#endif
