#ifndef TINKER_GPU_E_POLAR_H_
#define TINKER_GPU_E_POLAR_H_

#include "elec.h"
#include "energy_buffer.h"
#include "polgrp.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN elec_t epolar_electyp;

TINKER_EXTERN real u1scale, u2scale, u3scale, u4scale;
TINKER_EXTERN real d1scale, d2scale, d3scale, d4scale;
TINKER_EXTERN real p2scale, p3scale, p4scale, p5scale;
TINKER_EXTERN real p2iscale, p3iscale, p4iscale, p5iscale;

TINKER_EXTERN real udiag;

TINKER_EXTERN real* polarity;
TINKER_EXTERN real* thole;
TINKER_EXTERN real* pdamp;
TINKER_EXTERN real* polarity_inv;

TINKER_EXTERN NonbondedEnergy ep_handle;

TINKER_EXTERN DeviceVector<real, 3> ufld_vec;
TINKER_EXTERN DeviceVector<real, 6> dufld_vec;

TINKER_EXTERN DeviceVector<real, 3> work01_, work02_, work03_, work04_, work05_,
    work06_, work07_, work08_, work09_, work10_;

void epolar_data(rc_op op);

// see also subroutine epolar0e in epolar.f
void epolar0_dotprod(const real (*gpu_uind)[3], const real (*gpu_udirp)[3]);

// electrostatic field due to permanent multipoles
void dfield_coulomb(real* gpu_field, real* gpu_fieldp);
void dfield_ewald(real* gpu_field, real* gpu_fieldp);
void dfield_ewald_recip_self(real* gpu_field);
void dfield_ewald_real(real* gpu_field, real* gpu_fieldp);

// mutual electrostatic field due to induced dipole moments
void ufield_coulomb(const real* gpu_uind, const real* gpu_uinp, real* gpu_field,
                    real* gpu_fieldp);
void ufield_ewald(const real* gpu_uind, const real* gpu_uinp, real* gpu_field,
                  real* gpu_fieldp);
void ufield_ewald_recip_self(const real* gpu_uind, const real* gpu_uinp,
                             real* gpu_field, real* gpu_fieldp);
void ufield_ewald_real(const real* gpu_uind, const real* gpu_uinp,
                       real* gpu_field, real* gpu_fieldp);

void dfield(real* gpu_field, real* gpu_fieldp);
// -Tu operator
void ufield(const real* gpu_uind, const real* gpu_uinp, real* gpu_field,
            real* gpu_fieldp);

// different induction algorithms
void induce_mutual_pcg1(real* gpu_ud, real* gpu_up);
void induce(real* gpu_ud, real* gpu_up);

void epolar_coulomb(int vers);
void epolar_ewald(int vers);
void epolar(int vers);
TINKER_NAMESPACE_END

#endif
