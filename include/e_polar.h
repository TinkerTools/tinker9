#ifndef TINKER_GPU_E_POLAR_H_
#define TINKER_GPU_E_POLAR_H_

#include "dev_array.h"
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

TINKER_EXTERN device_pointer<real> polarity, thole, pdamp, polarity_inv;

TINKER_EXTERN NonbondedEnergy ep_handle;

TINKER_EXTERN device_pointer<real, 3> ufld;
TINKER_EXTERN device_pointer<real, 6> dufld;

TINKER_EXTERN device_pointer<real, 3> work01_, work02_, work03_, work04_,
    work05_, work06_, work07_, work08_, work09_, work10_;

void epolar_data(rc_op op);

// see also subroutine epolar0e in epolar.f
void epolar0_dotprod(const real (*uind)[3], const real (*udirp)[3]);

// electrostatic field due to permanent multipoles
void dfield_coulomb(real (*field)[3], real (*fieldp)[3]);
void dfield_ewald(real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_recip_self(real (*field)[3]);
void dfield_ewald_real(real (*field)[3], real (*fieldp)[3]);

// mutual electrostatic field due to induced dipole moments
void ufield_coulomb(const real (*uind)[3], const real (*uinp)[3],
                    real (*field)[3], real (*fieldp)[3]);
void ufield_ewald(const real (*uind)[3], const real (*uinp)[3],
                  real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_recip_self(const real (*uind)[3], const real (*uinp)[3],
                             real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real(const real (*uind)[3], const real (*uinp)[3],
                       real (*field)[3], real (*fieldp)[3]);

void dfield(real (*field)[3], real (*fieldp)[3]);
// -Tu operator
void ufield(const real (*uind)[3], const real (*uinp)[3], real (*field)[3],
            real (*fieldp)[3]);

// different induction algorithms
void induce_mutual_pcg1(real (*uind)[3], real (*uinp)[3]);
void induce(real (*uind)[3], real (*uinp)[3]);

void epolar_coulomb(int vers);
void epolar_ewald(int vers);
void epolar(int vers);
TINKER_NAMESPACE_END

#endif
