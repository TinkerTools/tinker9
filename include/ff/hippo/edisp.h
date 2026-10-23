#pragma once
#include "ff/energybuffer.h"
#include "ff/evdw.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup disp
/// \{
bool useDEwald();
void edispData(RcOp);
void edisp(int vers);

/// How the reciprocal and self terms of Ewald dispersion weight a mutated
/// ligand by its soft-core weight (edisp.f:dspmutwt). With \c split, the grid
/// holds the environment C6 in its real part and the ligand C6, the sites
/// \c part flags, in its imaginary part; the convolution separates the two
/// transforms from the copy \c qfwd of the forward transform.
struct DispSplit
{
   bool split;
   real wel; ///< ligand-environment weight
   real wll; ///< ligand-ligand weight
   const int* part;
   const real (*qfwd)[2];
};
/// \}
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
TINKER_EXTERN real csixpr;
TINKER_EXTERN real* csix;
TINKER_EXTERN real* adisp;

TINKER_EXTERN int ndspexclude;
TINKER_EXTERN int (*dspexclude)[2];
TINKER_EXTERN real* dspexclude_scale;

TINKER_EXTERN CountBuffer ndisp;
TINKER_EXTERN EnergyBuffer edsp;
TINKER_EXTERN VirialBuffer vir_edsp;
TINKER_EXTERN grad_prec* dedspx;
TINKER_EXTERN grad_prec* dedspy;
TINKER_EXTERN grad_prec* dedspz;
TINKER_EXTERN energy_prec energy_edsp;
TINKER_EXTERN virial_prec virial_edsp[9];

TINKER_EXTERN real dsp2scale;
TINKER_EXTERN real dsp3scale;
TINKER_EXTERN real dsp4scale;
TINKER_EXTERN real dsp5scale;

TINKER_EXTERN real elrc_vol_dsp;
TINKER_EXTERN real vlrc_vol_dsp;
}
