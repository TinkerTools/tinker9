#ifndef TINKER_PME_H_
#define TINKER_PME_H_

#include "gen_unit.h"
#include "rc_man.h"
#include <cmath>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * particle mesh ewald girds and parameters
 *
 * @code{.f}
 * !! allocate igrid(3,n)
 * !! allocate (bsbuild(bsorder,bsorder))
 * !! allocate (thetai1(4,bsorder,n))
 * !! allocate (thetai2(4,bsorder,n))
 * !! allocate (thetai3(4,bsorder,n))
 * !! allocate (qfac(nfft1,nfft2,nfft3))
 * allocate (bsmod1(nfft1))
 * allocate (bsmod2(nfft2))
 * allocate (bsmod3(nfft3))
 * allocate (qgrid(2,nfft1,nfft2,nfft3))
 * @endcode
 */
struct PME {
  struct Params {
    real aewald;
    int nfft1, nfft2, nfft3, bsorder;
    Params(real a, int n1, int n2, int n3, int o)
        : aewald(a)
        , nfft1(n1)
        , nfft2(n2)
        , nfft3(n3)
        , bsorder(o) {}
    bool operator==(const Params& st) const {
      const double eps = 1.0e-6;
      bool ans = std::fabs(aewald - st.aewald) < eps && nfft1 == st.nfft1 &&
          nfft2 == st.nfft2 && nfft3 == st.nfft3 && bsorder == st.bsorder;
      return ans;
    }
  };

  real aewald;
  int nfft1, nfft2, nfft3, bsorder;
  int* igrid;                     // deviceptr
  real *bsmod1, *bsmod2, *bsmod3; // deviceptr
  real* qgrid;                    // deviceptr

  bool operator==(const Params& p) const {
    Params p0(aewald, nfft1, nfft2, nfft3, bsorder);
    return p0 == p;
  }

  void set_params(const Params& p) {
    aewald = p.aewald;
    nfft1 = p.nfft1;
    nfft2 = p.nfft2;
    nfft3 = p.nfft3;
    bsorder = p.bsorder;
  }
};

typedef GenericUnit<PME, 1> PMEUnit;
TINKER_EXTERN PMEUnit epme_unit;  // electrostatic
TINKER_EXTERN PMEUnit ppme_unit;  // polarization
TINKER_EXTERN PMEUnit dpme_unit;  // dispersion
TINKER_EXTERN PMEUnit pvpme_unit; // polarization virial

TINKER_EXTERN double ewald_switch_cut, ewald_switch_off;

TINKER_EXTERN real (*cmp)[10];
TINKER_EXTERN real (*fmp)[10];
TINKER_EXTERN real (*cphi)[10];
TINKER_EXTERN real (*fphi)[20];

TINKER_EXTERN real (*fuind)[3];
TINKER_EXTERN real (*fuinp)[3];
TINKER_EXTERN real (*fdip_phi1)[10];
TINKER_EXTERN real (*fdip_phi2)[10];
TINKER_EXTERN real (*cphidp)[10];
TINKER_EXTERN real (*fphidp)[20];

TINKER_EXTERN real* vir_m;

void pme_data(rc_op op);

/// This function must be called after pme_data has been called because it
/// needs to know the number of pme objects created.
void fft_data(rc_op op);
void fftfront(PMEUnit pme_u);
void fftback(PMEUnit pme_u);

int use_ewald();
void pme_init(int vers);
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * make the scalar summation over reciprocal lattice
 */
void pme_conv0(PMEUnit pme_u);                 // without virial
void pme_conv1(PMEUnit pme_u, real* gpu_vir9); // with virial

void rpole_to_cmp();
/**
 * @brief
 * Input: cmp, cartesian rotated mpole.
 * Output: fmp, fractional rotated mpole.
 */
void cmp_to_fmp(PMEUnit pme_u, const real (*cmp)[10], real (*fmp)[10]);
void cuind_to_fuind(PMEUnit pme_u, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3]);
/**
 * @brief
 * Input: fphi.
 * Output: cphi.
 */
void fphi_to_cphi(PMEUnit pme_u, const real (*fphi)[20], real (*cphi)[10]);
/**
 * @brief
 * Input: fmp.
 * Output: qgrid.
 */
void grid_mpole(PMEUnit pme_u, real (*gpu_fmp)[10]);
/**
 * @brief
 * Input: qgrid.
 * Output: fphi.
 */
void fphi_mpole(PMEUnit pme_u, real (*gpu_fphi)[20]);

void grid_uind(PMEUnit pme_u, real (*gpu_find)[3], real (*gpu_finp)[3]);
void fphi_uind(PMEUnit pme_u, real (*gpu_fdip_phi1)[10],
               real (*gpu_fdip_phi2)[10], real (*gpu_fdip_sum_phi)[20]);
void fphi_uind2(PMEUnit pme_u, real (*gpu_fdip_phi1)[10],
                real (*gpu_fdip_phi2)[10]);
TINKER_NAMESPACE_END

#endif
