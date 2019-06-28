#include "gpu/decl_pme.h"
#include "gpu/decl_potent.h"
#include "gpu/decl_switch.h"
#include "gpu/e_mpole.h"
#include "gpu/rc.h"
#include "util/fort_str.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int empole_electyp;
std::string empole_electyp_str;

real* em;
int* nem;
real* vir_em;

void get_empole_type(int& typ, std::string& typ_str) {
  if (use_ewald()) {
    typ = elec_ewald;
    typ_str = "EWALD";
  } else {
    typ = elec_coulomb;
    typ_str = "COULOMB";
  }
}

void empole_data(rc_t rc) {
  if (!use_potent(mpole_term))
    return;

  if (rc & rc_dealloc)
    free_nev(nem, em, vir_em);

  if (rc & rc_alloc)
    alloc_nev(&nem, &em, &vir_em);

  if (rc & rc_copyin) {
    get_empole_type(empole_electyp, empole_electyp_str);

    if (empole_electyp == elec_coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);
  }
}

void empole(int vers) {
  if (empole_electyp == elec_coulomb)
    empole_coulomb(vers);
  else if (empole_electyp == elec_ewald)
    empole_ewald(vers);
}
}
TINKER_NAMESPACE_END
