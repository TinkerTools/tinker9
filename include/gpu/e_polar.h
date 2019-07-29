#ifndef TINKER_GPU_E_POLAR_H_
#define TINKER_GPU_E_POLAR_H_

#include "gpu/decl_elec.h"
#include "mod_polgrp.h"
#include "util_cxx.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern int epolar_electyp;
extern std::string epolar_electyp_str;

extern real u1scale, u2scale, u3scale, u4scale;
extern real d1scale, d2scale, d3scale, d4scale;
extern real p2scale, p3scale, p4scale, p5scale;
extern real p2iscale, p3iscale, p4iscale, p5iscale;

extern real* polarity;
extern real* thole;
extern real* pdamp;
extern real* polarity_inv;

extern real* ep;
extern int* nep;
extern real* vir_ep;

extern real (*ufld)[3];
extern real (*dufld)[6];

extern real (*work01_)[3];
extern real (*work02_)[3];
extern real (*work03_)[3];
extern real (*work04_)[3];
extern real (*work05_)[3];
extern real (*work06_)[3];
extern real (*work07_)[3];
extern real (*work08_)[3];
extern real (*work09_)[3];
extern real (*work10_)[3];

void get_epolar_type(int& typ, std::string& typ_str);
void epolar_data(rc_t rc);

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
}
TINKER_NAMESPACE_END

#endif
