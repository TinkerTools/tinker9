#ifndef TINKER_GPU_E_MPOLE_H_
#define TINKER_GPU_E_MPOLE_H_

#include "decl_elec.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern int empole_electyp;
extern std::string empole_electyp_str;

extern real* em;
extern int* nem;
extern real* vir_em;

void get_empole_type(int& typ, std::string& typ_str);
void empole_data(rc_t rc);

void empole_coulomb(int vers);
void empole_ewald(int vers);
void empole(int vers);
}
TINKER_NAMESPACE_END

#endif
