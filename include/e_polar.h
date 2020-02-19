#pragma once
#include "dev_array.h"
#include "elec.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN elec_t epolar_electyp;

TINKER_EXTERN real u1scale, u2scale, u3scale, u4scale;
TINKER_EXTERN real d1scale, d2scale, d3scale, d4scale;
TINKER_EXTERN real p2scale, p3scale, p4scale, p5scale;
TINKER_EXTERN real p2iscale, p3iscale, p4iscale, p5iscale;

TINKER_EXTERN int nuexclude_;
TINKER_EXTERN device_pointer<int, 2> uexclude_;
TINKER_EXTERN device_pointer<real> uexclude_scale_;

TINKER_EXTERN int ndpexclude_;
TINKER_EXTERN device_pointer<int, 2> dpexclude_;
TINKER_EXTERN device_pointer<real, 2> dpexclude_scale_;

TINKER_EXTERN int ndpuexclude_;
TINKER_EXTERN device_pointer<int, 2> dpuexclude_;
TINKER_EXTERN device_pointer<real, 3> dpuexclude_scale_;

TINKER_EXTERN real udiag;

TINKER_EXTERN device_pointer<real> polarity, thole, pdamp, polarity_inv;

TINKER_EXTERN count_buffer nep;
TINKER_EXTERN energy_buffer ep;
TINKER_EXTERN virial_buffer vir_ep;

TINKER_EXTERN device_pointer<real, 3> ufld;
TINKER_EXTERN device_pointer<real, 6> dufld;

TINKER_EXTERN device_pointer<real, 3> work01_, work02_, work03_, work04_,
   work05_, work06_, work07_, work08_, work09_, work10_;

void epolar_data(rc_op op);


// different induction algorithms
void induce_mutual_pcg1(real (*uind)[3], real (*uinp)[3]);
void induce_mutual_pcg1_cu(real (*uind)[3], real (*uinp)[3]);
void induce(real (*uind)[3], real (*uinp)[3]);


// electrostatic field due to permanent multipoles
// clang-format off
void dfield                  (real (*field)[3], real (*fieldp)[3]);
void dfield_nonewald         (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald            (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_recip_self (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_real       (real (*field)[3], real (*fieldp)[3]);

void dfield_nonewald_acc         (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_recip_self_acc (real (*field)[3]);
void dfield_ewald_real_acc       (real (*field)[3], real (*fieldp)[3]);
void dfield_nonewald_cu          (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_real_cu        (real (*field)[3], real (*fieldp)[3]);
// clang-format on


// mutual electrostatic field due to induced dipole moments
// clang-format off
// -Tu operator
void ufield                  (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_nonewald         (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald            (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_recip_self (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real       (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);

void ufield_nonewald_acc         (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_recip_self_acc (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real_acc       (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_nonewald_cu          (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real_cu        (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
// clang-format on


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
