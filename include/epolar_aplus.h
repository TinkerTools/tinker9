#pragma once
#include "epolar.h"
#include "mod.polar.h"


namespace tinker {
void epolar_aplus_data(rc_op op);


// different induction algorithms
void induce_mutual_pcg3(real (*uind)[3]);
void induce_mutual_pcg3_acc(real (*uind)[3]);
void induce_mutual_pcg3_cu(real (*uind)[3]);
void induce3(real (*uind)[3]);


void epolar_aplus(int vers);
void epolar_aplus_nonewald(int vers, int use_cf);
void epolar_aplus_ewald(int vers, int use_cf);
void epolar_aplus_ewald_real(int vers, int use_cf);
void epolar_aplus_ewald_recip_self(int vers, int use_cf);

// see also subroutine epolar0e in epolar.f

void epolar_aplus_nonewald_acc(int vers, int use_cf, const real (*d)[3]);
void epolar_aplus_ewald_real_acc(int vers, int use_cf, const real (*d)[3]);
void epolar_aplus_ewald_recip_self_acc(int vers, int use_cf,
                                        const real (*d)[3]);
void epolar_aplus_nonewald_cu(int vers, int use_cf, const real (*d)[3]);
void epolar_aplus_ewald_real_cu(int vers, int use_cf, const real (*d)[3]);
}
