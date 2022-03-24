#pragma once
#include "ff/elec.h"
#include "ff/pmestuf.h"
#include "mod/elecamoeba.h"

namespace tinker {
void epolar_data(RcOp);

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
void epolar_ewald_recip_self_acc(int vers, const real (*d)[3], const real (*p)[3]);
void epolar0_dotprod_acc(const real (*uind)[3], const real (*udirp)[3]);
void epolar_nonewald_cu(int vers, const real (*d)[3], const real (*p)[3]);
void epolar_ewald_real_cu(int vers, const real (*d)[3], const real (*p)[3]);
}
