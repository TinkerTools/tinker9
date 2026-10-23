#pragma once
#include "ff/amoeba/mpole.h"
#include "ff/energybuffer.h"
#include "ff/evdw.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup repel
/// Soft core scaling of a mutated Pauli repulsion pair,
/// \f$ s = \lambda^n r / \sqrt{\alpha \lambda^m (1-\lambda)^p + r^2} \f$,
/// as the scale \c rlsc and shift \c rshf, with the lambda derivative factors
/// of erepsc.f:repsc4 (\c rshf1 and \c rshf2 carry lambda and lambda squared
/// times the shift derivatives) and the van der Waals lambda chain rule.
struct RepelSoft
{
   real rlsc, rshf, n;
   real rlsc1, rlsc2;
   real rshf1, rshf2;
   real dvl, d2vl;
};

/// \ingroup repel
/// The soft core terms at the current vlam (erepsc.f:repsc). The lambda
/// derivative factors are set only if \c dodl, and the second derivative ones
/// only with use_d2lmda (erepsc.f:repsc4).
RepelSoft repelSoft(bool dodl);

/// \ingroup repel
void erepelData(RcOp);
/// \ingroup repel
void erepel(int vers);
/// \ingroup repel
void repoleInit(int vers);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
TINKER_EXTERN real (*repole)[MPL_TOTAL];
TINKER_EXTERN real (*rrepole)[MPL_TOTAL];
TINKER_EXTERN real* sizpr;
TINKER_EXTERN real* dmppr;
TINKER_EXTERN real* elepr;

TINKER_EXTERN int nrepexclude;
TINKER_EXTERN int (*repexclude)[2];
TINKER_EXTERN real* repexclude_scale;

TINKER_EXTERN CountBuffer nrep;
TINKER_EXTERN EnergyBuffer er;
TINKER_EXTERN VirialBuffer vir_er;
TINKER_EXTERN grad_prec* derx;
TINKER_EXTERN grad_prec* dery;
TINKER_EXTERN grad_prec* derz;
TINKER_EXTERN energy_prec energy_er;
TINKER_EXTERN virial_prec virial_er[9];

TINKER_EXTERN real r2scale;
TINKER_EXTERN real r3scale;
TINKER_EXTERN real r4scale;
TINKER_EXTERN real r5scale;

TINKER_EXTERN double scalphar;
TINKER_EXTERN double scexpr;
TINKER_EXTERN double scexprm;
TINKER_EXTERN double scexprp;
}
