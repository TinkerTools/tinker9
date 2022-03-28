#pragma once
#include "tool/energybuffer.h"
#include "tool/fft.h"
#include "tool/genunit.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup pme
/// \brief Particle mesh ewald girds and parameters.
///
/// \code{.f}
/// !! In Tinker:
/// allocate igrid(3,n)
/// allocate (bsbuild(bsorder,bsorder))
/// allocate (thetai1(4,bsorder,n))
/// allocate (thetai2(4,bsorder,n))
/// allocate (thetai3(4,bsorder,n))
/// allocate (qfac(nfft1,nfft2,nfft3))
/// allocate (bsmod1(nfft1))
/// allocate (bsmod2(nfft2))
/// allocate (bsmod3(nfft3))
/// allocate (qgrid(2,nfft1,nfft2,nfft3))
/// \endcode
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
using PMEUnit = GenericUnit<PME, GenericUnitVersion::ENABLE_ON_DEVICE>;
}

namespace tinker {
void pmeData(RcOp);
// This function must be called after pmeData has been called because it
// needs to know the number of pme objects created.
void fftData(RcOp);
void fftfront(PMEUnit pme_u);
void fftback(PMEUnit pme_u);
}

namespace tinker {
void bspline_fill(PMEUnit, int level);

void grid_pchg(PMEUnit, real* pchg);
void grid_disp(PMEUnit, real* csix);
void grid_mpole(PMEUnit, real (*fmp)[10]);
void grid_uind(PMEUnit, real (*find)[3], real (*finp)[3]);

void pme_conv(PMEUnit);                 // update grid
void pme_conv(PMEUnit, VirialBuffer v); // update grid and accumulate vterm
void pme_conv(PMEUnit, EnergyBuffer e); // update grid and accumulate eterm
void pme_conv(PMEUnit, EnergyBuffer e, VirialBuffer v);

void fphi_mpole(PMEUnit);
void fphi_uind(PMEUnit, real (*fdip_phi1)[10], real (*fdip_phi2)[10], real (*fdip_sum_phi)[20]);
void fphi_uind2(PMEUnit, real (*fdip_phi1)[10], real (*fdip_phi2)[10]);

void rpole_to_cmp();
void cmp_to_fmp(PMEUnit, const real (*cmp)[10], real (*fmp)[10]);
void cuind_to_fuind(
   PMEUnit, const real (*cind)[3], const real (*cinp)[3], real (*fuind)[3], real (*fuinp)[3]);
void fphi_to_cphi(PMEUnit, const real (*fphi)[20], real (*cphi)[10]);

//====================================================================//

void bspline_fill_cu(PMEUnit, int level);

#define TINKER_CU_THETA_ON_THE_FLY_GRID_MPOLE 1
#define TINKER_CU_THETA_ON_THE_FLY_GRID_UIND  0

void grid_pchg_acc(PMEUnit, real*);
void grid_pchg_cu(PMEUnit, real*);
void grid_disp_acc(PMEUnit, real*);
void grid_disp_cu(PMEUnit, real*);
void grid_mpole_acc(PMEUnit, real (*)[10]);
void grid_mpole_cu(PMEUnit, real (*)[10]);
void grid_uind_acc(PMEUnit, real (*)[3], real (*)[3]);
void grid_uind_cu(PMEUnit, real (*)[3], real (*)[3]);

void pme_conv_acc(PMEUnit, EnergyBuffer, VirialBuffer);
void pme_conv_cu(PMEUnit, EnergyBuffer, VirialBuffer);

void fphi_mpole_acc(PMEUnit, real (*)[20]);
void fphi_mpole_cu(PMEUnit, real (*)[20]);
void fphi_uind_acc(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
void fphi_uind_cu(PMEUnit, real (*)[10], real (*)[10], real (*)[20]);
void fphi_uind2_acc(PMEUnit, real (*)[10], real (*)[10]);
void fphi_uind2_cu(PMEUnit, real (*)[10], real (*)[10]);

void rpole_to_cmp_acc();
void cmp_to_fmp_acc(PMEUnit, const real (*)[10], real (*)[10]);
void cuind_to_fuind_acc(PMEUnit, const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
void fphi_to_cphi_acc(PMEUnit, const real (*)[20], real (*)[10]);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
TINKER_EXTERN PMEUnit epme_unit;  // electrostatic
TINKER_EXTERN PMEUnit ppme_unit;  // polarization
TINKER_EXTERN PMEUnit pvpme_unit; // polarization virial
TINKER_EXTERN PMEUnit dpme_unit;  // dispersion

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

TINKER_EXTERN VirialBuffer vir_m;
}
