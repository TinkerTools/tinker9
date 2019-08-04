#ifndef TINKER_GPU_E_MPOLE_H_
#define TINKER_GPU_E_MPOLE_H_

#include "cxx.h"
#include "mod_elec.h"

TINKER_NAMESPACE_BEGIN
extern int empole_electyp;
extern std::string empole_electyp_str;

extern real m2scale, m3scale, m4scale, m5scale;

extern real* em;
extern int* nem;
extern real* vir_em;

void get_empole_type(int& typ, std::string& typ_str);
void empole_data(rc_op op);

void empole_coulomb(int vers);
void empole_ewald(int vers);
void empole(int vers);
TINKER_NAMESPACE_END

#endif
