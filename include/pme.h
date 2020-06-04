#pragma once
#include "energy_buffer.h"
#include "fft.h"
#include "tool/darray.h"
#include "tool/gen_unit.h"
#include "tool/rc_man.h"


namespace tinker {
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


void pme_data(rc_op op);
// This function must be called after pme_data has been called because it
// needs to know the number of pme objects created.
void fft_data(rc_op op);
void fftfront(PMEUnit pme_u);
void fftback(PMEUnit pme_u);
using FFTPlanUnit = GenericUnit<FFTPlan, GenericUnitVersion::DisableOnDevice>;
}
