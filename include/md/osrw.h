#pragma once
#include "ff/amoeba/mpole.h"
#include "ff/precision.h"
#include "ff/timescale.h"
#include "tool/rcman.h"

namespace tinker {
void osrw_mech();
void osrwData(RcOp);
double osrw_lam_expr0(int form, double lam);
double osrw_lam_expr1(int form, double lam);
double osrw_lam_expr2(int form, double lam);
void osrw_altele(double);
void osrw_alttor(double);
void osrw_altvdw(double);
void osrw_energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);
void osrw_energy(int vers);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
constexpr int OSRW_LAM_LINEAR = 0;
constexpr int OSRW_LAM_QUADRATIC = 1;

TINKER_EXTERN bool use_osrw;
TINKER_EXTERN double osrw_lambda;
TINKER_EXTERN int osrw_vdw;
TINKER_EXTERN int osrw_ele;
TINKER_EXTERN int osrw_tor;

/**
 * \ingroup osrw
 * \{
 * \var osrw_du1
 * \brief Partial derivative of potential energy w.r.t. lambda.
 * \var osrw_dv1
 * \brief Partial derivative of virial tensor w.r.t. lambda.
 * \var osrw_dgx
 * \brief Partial derivative of energy gradient w.r.t. lambda.
 * \var osrw_dgy
 * \copydoc osrw_dgx
 * \var osrw_dgz
 * \copydoc osrw_dgx
 * \}
 */
TINKER_EXTERN energy_prec osrw_du1;
TINKER_EXTERN virial_prec osrw_dv1[9];
TINKER_EXTERN grad_prec *osrw_dgx, *osrw_dgy, *osrw_dgz;

TINKER_EXTERN real* osrw_pchg;
TINKER_EXTERN real (*osrw_pole)[MPL_TOTAL];
TINKER_EXTERN real* osrw_polarity;
TINKER_EXTERN int osrw_ntbnd;
TINKER_EXTERN int (*osrw_itbnd)[2];
TINKER_EXTERN real (*osrw_tors1)[4];
TINKER_EXTERN real (*osrw_tors2)[4];
TINKER_EXTERN real (*osrw_tors3)[4];
TINKER_EXTERN real (*osrw_tors4)[4];
TINKER_EXTERN real (*osrw_tors5)[4];
TINKER_EXTERN real (*osrw_tors6)[4];
}
