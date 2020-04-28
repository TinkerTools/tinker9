#pragma once
#include "elec.h"
#include "energy_buffer.h"
#include "pmestuf.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN real u1scale, u2scale, u3scale, u4scale;
TINKER_EXTERN real d1scale, d2scale, d3scale, d4scale;
TINKER_EXTERN real p2scale, p3scale, p4scale, p5scale;
TINKER_EXTERN real p2iscale, p3iscale, p4iscale, p5iscale;

TINKER_EXTERN int nuexclude;
TINKER_EXTERN pointer<int, 2> uexclude;
TINKER_EXTERN pointer<real> uexclude_scale;

TINKER_EXTERN int ndpexclude;
TINKER_EXTERN pointer<int, 2> dpexclude;
TINKER_EXTERN pointer<real, 2> dpexclude_scale;

TINKER_EXTERN int ndpuexclude;
TINKER_EXTERN pointer<int, 2> dpuexclude;
TINKER_EXTERN pointer<real, 3> dpuexclude_scale;

TINKER_EXTERN real udiag;

TINKER_EXTERN pointer<real> polarity, thole, pdamp, polarity_inv;

TINKER_EXTERN count_buffer nep;
TINKER_EXTERN energy_buffer ep;
TINKER_EXTERN virial_buffer vir_ep;
TINKER_EXTERN grad_prec *depx, *depy, *depz;

TINKER_EXTERN pointer<real, 3> ufld;
TINKER_EXTERN pointer<real, 6> dufld;

TINKER_EXTERN pointer<real, 3> work01_, work02_, work03_, work04_, work05_,
   work06_, work07_, work08_, work09_, work10_;

void epolar_data(rc_op op);


// different induction algorithms
void induce_mutual_pcg1(real (*uind)[3], real (*uinp)[3]);
void induce_mutual_pcg1_acc(real (*uind)[3], real (*uinp)[3]);
void induce_mutual_pcg1_cu(real (*uind)[3], real (*uinp)[3]);
void induce(real (*uind)[3], real (*uinp)[3]);


void epolar(int vers);
void epolar_nonewald(int vers);
void epolar_ewald(int vers);
void epolar_ewald_real(int vers);
void epolar_ewald_recip_self(int vers);
// see also subroutine epolar0e in epolar.f
void epolar0_dotprod(const real (*uind)[3], const real (*udirp)[3]);

void epolar_nonewald_acc(int vers, const real (*d)[3], const real (*p)[3]);
void epolar_ewald_real_acc(int vers, const real (*d)[3], const real (*p)[3]);
void epolar_ewald_recip_self_acc(int vers, const real (*d)[3],
                                 const real (*p)[3]);
void epolar0_dotprod_acc(const real (*uind)[3], const real (*udirp)[3]);
void epolar_nonewald_cu(int vers, const real (*d)[3], const real (*p)[3]);
void epolar_ewald_real_cu(int vers, const real (*d)[3], const real (*p)[3]);
TINKER_NAMESPACE_END
