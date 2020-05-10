#pragma once
#include "elec.h"
#include "energy.h"
#include "rc_man.h"


namespace tinker {
constexpr int OSRW_LAM_LINEAR = 0;
constexpr int OSRW_LAM_QUADRATIC = 1;


extern bool use_osrw;
extern double osrw_lambda;
extern int osrw_vdw;
extern int osrw_ele;
extern int osrw_tor;


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
extern energy_prec osrw_du1;
extern virial_prec osrw_dv1[9];
extern grad_prec *osrw_dgx, *osrw_dgy, *osrw_dgz;


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
void osrw_altele_acc(double);
void osrw_alttor_acc(double);
void osrw_energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);
void osrw_energy(int vers);
}
