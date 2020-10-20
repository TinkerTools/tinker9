#pragma once
#include "elec.h"
#include "mod.chgpen.h"
#include "mod.polar.h"
#include "mod.polpot.h"
#include "pmestuf.h"


namespace tinker {
void epolar_chgpen_data(rc_op op);


// different induction algorithms
void induce_mutual_pcg2(real (*uind)[3]);
void induce_mutual_pcg_acc2(real (*uind)[3]);
void induce_mutual_pcg_cu2(real (*uind)[3]);
void induce2(real (*uind)[3]);


void epolar_chgpen(int vers);
void epolar_chgpen_nonewald(int vers, int use_cf);
void epolar_chgpen_ewald(int vers, int use_cf);
void epolar_chgpen_ewald_real(int vers, int use_cf);
void epolar_chgpen_ewald_recip_self(int vers, int use_cf);
void epolar1_dotprod(const real (*uind)[3], const real (*udirp)[3]);

// see also subroutine epolar0e in epolar.f

void epolar1_dotprod_acc(const real (*uind)[3], const real (*udirp)[3]);
void epolar_chgpen_nonewald_acc(int vers, const real (*d)[3]);
void epolar_chgpen_ewald_real_acc(int vers, const real (*d)[3]);
void epolar_chgpen_ewald_recip_self_acc(int vers, int use_cf, const real (*d)[3]);
void epolar_chgpen_nonewald_cu(int vers, int use_cf, const real (*d)[3]);
void epolar_chgpen_ewald_real_cu(int vers, int use_cf, const real (*d)[3]);
}
