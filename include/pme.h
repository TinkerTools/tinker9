#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "fft.h"
#include "gen_unit.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup pme
 * \brief Particle mesh ewald girds and parameters.
 *
 * \code{.f}
 * !! In Tinker:
 * allocate igrid(3,n)
 * allocate (bsbuild(bsorder,bsorder))
 * allocate (thetai1(4,bsorder,n))
 * allocate (thetai2(4,bsorder,n))
 * allocate (thetai3(4,bsorder,n))
 * allocate (qfac(nfft1,nfft2,nfft3))
 * allocate (bsmod1(nfft1))
 * allocate (bsmod2(nfft2))
 * allocate (bsmod3(nfft3))
 * allocate (qgrid(2,nfft1,nfft2,nfft3))
 * \endcode
 */
struct PME
{
   real aewald;
   int nfft1, nfft2, nfft3, bsorder;
   real *bsmod1, *bsmod2, *bsmod3;
   real* qgrid;
   int* igrid;
   real *thetai1, *thetai2, *thetai3;

   struct Params
   {
      real aewald;
      int nfft1, nfft2, nfft3, bsorder;
      bool operator==(const Params& st) const;
      Params(real a, int n1, int n2, int n3, int o);
   };
   void set_params(const Params& p);
   PME::Params get_params() const;
   bool operator==(const Params& p) const;

   ~PME();
};
using PMEUnit = GenericUnit<PME, GenericUnitVersion::EnableOnDevice>;
TINKER_EXTERN PMEUnit epme_unit;  // electrostatic
TINKER_EXTERN PMEUnit ppme_unit;  // polarization
TINKER_EXTERN PMEUnit dpme_unit;  // dispersion
TINKER_EXTERN PMEUnit pvpme_unit; // polarization virial

TINKER_EXTERN pointer<real, 10> cmp, fmp, cphi;
TINKER_EXTERN pointer<real, 20> fphi;

TINKER_EXTERN pointer<real, 3> fuind, fuinp;
TINKER_EXTERN pointer<real, 10> fdip_phi1, fdip_phi2, cphidp;
TINKER_EXTERN pointer<real, 20> fphidp;

TINKER_EXTERN virial_buffer vir_m;

void pme_data(rc_op op);
// This function must be called after pme_data has been called because it
// needs to know the number of pme objects created.
void fft_data(rc_op op);

void fftfront(PMEUnit pme_u);
void fftback(PMEUnit pme_u);
using FFTPlanUnit = GenericUnit<FFTPlan, GenericUnitVersion::DisableOnDevice>;

void pme_init(int vers);
TINKER_NAMESPACE_END


TINKER_NAMESPACE_BEGIN
/// \ingroup pme
/// \brief Make the scalar summation over reciprocal lattice without virial.
void pme_conv(PMEUnit pme_u);
/// \ingroup pme
/// \brief Make the scalar summation over reciprocal lattice with virial.
void pme_conv(PMEUnit pme_u, virial_buffer v);
void pme_conv_acc(PMEUnit, virial_buffer);


/// \ingroup pme
void rpole_to_cmp();
/// \ingroup pme
void rpole_to_cmp_acc();


/**
 * \ingroup pme
 * \brief Get PME B-spline coefficients.
 *
 * Input: x, y, z, recip box
 * Output: igrid, thetai1, thetai2, thetai3
 */
void bspline_fill_cu(PMEUnit, int level);


/**
 * \ingroup pme
 *
 * Input: cmp, cartesian rotated mpole.
 * Output: fmp, fractional rotated mpole.
 */
void cmp_to_fmp(PMEUnit pme_u, const real (*cmp)[10], real (*fmp)[10]);
/**
 * \ingroup pme
 */
void cuind_to_fuind(PMEUnit pme_u, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3]);
/**
 * \ingroup pme
 * Input: fphi.
 * Output: cphi.
 */
void fphi_to_cphi(PMEUnit pme_u, const real (*fphi)[20], real (*cphi)[10]);
void cmp_to_fmp_acc(PMEUnit, const real (*)[10], real (*)[10]);
void cuind_to_fuind_acc(PMEUnit, const real (*)[3], const real (*)[3],
                        real (*)[3], real (*)[3]);
void fphi_to_cphi_acc(PMEUnit, const real (*)[20], real (*)[10]);


// grid_*: input: fmp; output: qgrid.
// fphi_*: input: qgrid; output: fphi.

/// \ingroup pme
void grid_mpole(PMEUnit pme_u, real (*gpu_fmp)[10]);
/// \ingroup pme
void grid_uind(PMEUnit pme_u, real (*gpu_find)[3], real (*gpu_finp)[3]);

void grid_mpole_acc(PMEUnit, real (*)[10]);
void grid_uind_acc(PMEUnit, real (*)[3], real (*)[3]);
void grid_mpole_cu(PMEUnit, real (*)[10]);
void grid_uind_cu(PMEUnit, real (*)[3], real (*)[3]);


/// \ingroup pme
void fphi_mpole(PMEUnit pme_u, real (*gpu_fphi)[20]);
/// \ingroup pme
void fphi_uind(PMEUnit pme_u, real (*gpu_fdip_phi1)[10],
               real (*gpu_fdip_phi2)[10], real (*gpu_fdip_sum_phi)[20]);
/// \ingroup pme
void fphi_uind2(PMEUnit pme_u, real (*gpu_fdip_phi1)[10],
                real (*gpu_fdip_phi2)[10]);

void fphi_mpole_acc(PMEUnit, real (*)[20]);
void fphi_uind_acc(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
void fphi_uind2_acc(PMEUnit pme_u, real (*)[10], real (*)[10]);
void fphi_mpole_cu(PMEUnit, real (*)[20]);
void fphi_uind_cu(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
void fphi_uind2_cu(PMEUnit pme_u, real (*)[10], real (*)[10]);


void pme_cuda_func_config();
TINKER_NAMESPACE_END
