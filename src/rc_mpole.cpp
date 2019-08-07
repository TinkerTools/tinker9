#include "gpu/e_mpole.h"
#include "io_fort_str.h"
#include "md.h"
#include "pme.h"
#include "potent.h"
#include "switch.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
int empole_electyp;
std::string empole_electyp_str;

real m2scale, m3scale, m4scale, m5scale;

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

void empole_data(rc_op op) {
  if (!use_potent(mpole_term))
    return;

  if (op & rc_dealloc)
    dealloc_nev(nem, em, vir_em);

  if (op & rc_alloc)
    alloc_nev(&nem, &em, &vir_em);

  if (op & rc_init) {
    get_empole_type(empole_electyp, empole_electyp_str);

    if (empole_electyp == elec_coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);

    m2scale = mplpot::m2scale;
    m3scale = mplpot::m3scale;
    m4scale = mplpot::m4scale;
    m5scale = mplpot::m5scale;
  }
}

void empole(int vers) {
  if (empole_electyp == elec_coulomb)
    empole_coulomb(vers);
  else if (empole_electyp == elec_ewald)
    empole_ewald(vers);
}
TINKER_NAMESPACE_END
