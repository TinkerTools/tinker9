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
TINKER_NAMESPACE_END
