#ifndef TINKER_MOD_POTENT_HH_
#define TINKER_MOD_POTENT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace potent {
extern int& use_bond;
extern int& use_angle;
extern int& use_strbnd;
extern int& use_urey;
extern int& use_angang;
extern int& use_opbend;
extern int& use_opdist;
extern int& use_improp;
extern int& use_imptor;
extern int& use_tors;
extern int& use_pitors;
extern int& use_strtor;
extern int& use_angtor;
extern int& use_tortor;
extern int& use_vdw;
extern int& use_repuls;
extern int& use_disp;
extern int& use_charge;
extern int& use_chgdpl;
extern int& use_dipole;
extern int& use_mpole;
extern int& use_polar;
extern int& use_chgtrn;
extern int& use_rxnfld;
extern int& use_solv;
extern int& use_metal;
extern int& use_geom;
extern int& use_extra;
extern int& use_born;
extern int& use_orbit;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(potent, use_bond);
extern "C" int m_tinker_mod(potent, use_angle);
extern "C" int m_tinker_mod(potent, use_strbnd);
extern "C" int m_tinker_mod(potent, use_urey);
extern "C" int m_tinker_mod(potent, use_angang);
extern "C" int m_tinker_mod(potent, use_opbend);
extern "C" int m_tinker_mod(potent, use_opdist);
extern "C" int m_tinker_mod(potent, use_improp);
extern "C" int m_tinker_mod(potent, use_imptor);
extern "C" int m_tinker_mod(potent, use_tors);
extern "C" int m_tinker_mod(potent, use_pitors);
extern "C" int m_tinker_mod(potent, use_strtor);
extern "C" int m_tinker_mod(potent, use_angtor);
extern "C" int m_tinker_mod(potent, use_tortor);
extern "C" int m_tinker_mod(potent, use_vdw);
extern "C" int m_tinker_mod(potent, use_repuls);
extern "C" int m_tinker_mod(potent, use_disp);
extern "C" int m_tinker_mod(potent, use_charge);
extern "C" int m_tinker_mod(potent, use_chgdpl);
extern "C" int m_tinker_mod(potent, use_dipole);
extern "C" int m_tinker_mod(potent, use_mpole);
extern "C" int m_tinker_mod(potent, use_polar);
extern "C" int m_tinker_mod(potent, use_chgtrn);
extern "C" int m_tinker_mod(potent, use_rxnfld);
extern "C" int m_tinker_mod(potent, use_solv);
extern "C" int m_tinker_mod(potent, use_metal);
extern "C" int m_tinker_mod(potent, use_geom);
extern "C" int m_tinker_mod(potent, use_extra);
extern "C" int m_tinker_mod(potent, use_born);
extern "C" int m_tinker_mod(potent, use_orbit);

int& use_bond = m_tinker_mod(potent, use_bond);
int& use_angle = m_tinker_mod(potent, use_angle);
int& use_strbnd = m_tinker_mod(potent, use_strbnd);
int& use_urey = m_tinker_mod(potent, use_urey);
int& use_angang = m_tinker_mod(potent, use_angang);
int& use_opbend = m_tinker_mod(potent, use_opbend);
int& use_opdist = m_tinker_mod(potent, use_opdist);
int& use_improp = m_tinker_mod(potent, use_improp);
int& use_imptor = m_tinker_mod(potent, use_imptor);
int& use_tors = m_tinker_mod(potent, use_tors);
int& use_pitors = m_tinker_mod(potent, use_pitors);
int& use_strtor = m_tinker_mod(potent, use_strtor);
int& use_angtor = m_tinker_mod(potent, use_angtor);
int& use_tortor = m_tinker_mod(potent, use_tortor);
int& use_vdw = m_tinker_mod(potent, use_vdw);
int& use_repuls = m_tinker_mod(potent, use_repuls);
int& use_disp = m_tinker_mod(potent, use_disp);
int& use_charge = m_tinker_mod(potent, use_charge);
int& use_chgdpl = m_tinker_mod(potent, use_chgdpl);
int& use_dipole = m_tinker_mod(potent, use_dipole);
int& use_mpole = m_tinker_mod(potent, use_mpole);
int& use_polar = m_tinker_mod(potent, use_polar);
int& use_chgtrn = m_tinker_mod(potent, use_chgtrn);
int& use_rxnfld = m_tinker_mod(potent, use_rxnfld);
int& use_solv = m_tinker_mod(potent, use_solv);
int& use_metal = m_tinker_mod(potent, use_metal);
int& use_geom = m_tinker_mod(potent, use_geom);
int& use_extra = m_tinker_mod(potent, use_extra);
int& use_born = m_tinker_mod(potent, use_born);
int& use_orbit = m_tinker_mod(potent, use_orbit);
#endif

} TINKER_NAMESPACE_END

#endif
