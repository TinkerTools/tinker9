#pragma once
#include "elec.h"
#include "energy.h"
#include "rc_man.h"


namespace tinker {
/**
 * \page kosrw  OSRW Keywords
 *
 * ### OSRW-ELE [LINEAR/QUADRATIC]
 *
 * ### OSRW-LAMBDA [real]
 * This keyword sets the internal logical flag for OSRW to `true` and provides
 * the initial value of lambda.
 *
 * ### OSRW-TORS [LINEAR/QUADRATIC]
 *
 * ### OSRW-VDW [LINEAR/QUADRATIC]
 *
 * ### ROTATABLE-BOND [integer list]
 */


constexpr int OSRW_LAM_LINEAR = 0;
constexpr int OSRW_LAM_QUADRATIC = 1;


extern bool use_osrw;
extern double osrw_lambda;
extern int osrw_vdw;
extern int osrw_ele;
extern int osrw_tor;


extern energy_prec osrw_du1;
extern virial_prec osrw_dv1[9];
extern grad_prec *osrw_dgx, *osrw_dgy, *osrw_dgz;
extern grad_prec *osrw_gx, *osrw_gy, *osrw_gz;
extern real* osrw_pchg;
extern real (*osrw_pole)[mpl_total];
extern real* osrw_polarity;
extern int osrw_ntbnd;
extern int (*osrw_itbnd)[2];
extern real (*osrw_tors1)[4];
extern real (*osrw_tors2)[4];
extern real (*osrw_tors3)[4];
extern real (*osrw_tors4)[4];
extern real (*osrw_tors5)[4];
extern real (*osrw_tors6)[4];


void osrw_mech();
void osrw_data(rc_op);
double osrw_lam_expr0(int form, double lam);
double osrw_lam_expr1(int form, double lam);
double osrw_lam_expr2(int form, double lam);
void osrw_altele(double);
void osrw_alttor(double);
void osrw_altvdw(double);
// U(L) = [1-h(L)] U0 + h(L) U1
// U' = -h'(L) U0 + h'(L) U1
void osrw_energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);
void osrw_energy(int vers);


void osrw_altele_acc(double);
void osrw_alttor_acc(double);
}
