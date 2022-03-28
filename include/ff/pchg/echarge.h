#pragma once
#include "tool/energybuffer.h"
#include "tool/rcman.h"

namespace tinker {

void echargeData(RcOp);

void echarge(int vers);

void echarge_nonewald(int);
void echarge_nonewald_acc(int);
void echarge_nonewald_cu(int);

void echarge_ewald_recip_self(int);
void echarge_ewald_fphi_self_acc(int);
void echarge_ewald_fphi_self_cu(int);

// void echarge_ewald_real(int);
void echarge_ewald_real_acc(int);
void echarge_ewald_real_cu(int);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

// charge chgpot
namespace tinker {
/// \ingroup charge
/// \brief Magnitude of the partial charges (e-) of each **atom**.
/// \note Unlike Tinker, where only non-zero charges will be stored, this
/// array also includes zero charges.
TINKER_EXTERN real* pchg;

TINKER_EXTERN CountBuffer nec;
TINKER_EXTERN EnergyBuffer ec;
TINKER_EXTERN VirialBuffer vir_ec;
TINKER_EXTERN grad_prec* decx;
TINKER_EXTERN grad_prec* decy;
TINKER_EXTERN grad_prec* decz;
TINKER_EXTERN energy_prec energy_ec;
TINKER_EXTERN virial_prec virial_ec[9];

TINKER_EXTERN real ebuffer;
TINKER_EXTERN real c2scale, c3scale, c4scale, c5scale;
TINKER_EXTERN int ncexclude;
TINKER_EXTERN int (*cexclude)[2];
TINKER_EXTERN real* cexclude_scale;
}
