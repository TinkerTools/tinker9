#include "e_mpole.h"
#include "ext/tinker/detail/mplpot.hh"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include "switch.h"

TINKER_NAMESPACE_BEGIN
void empole_data(rc_op op) {
  if (!use_potent(mpole_term))
    return;

  if (op & rc_dealloc)
    em_handle.dealloc();

  if (op & rc_alloc)
    em_handle.alloc(n);

  if (op & rc_init) {
    if (use_ewald()) {
      empole_electyp = elec_t::ewald;
    } else {
      empole_electyp = elec_t::coulomb;
    }

    if (empole_electyp == elec_t::coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);

    m2scale = mplpot::m2scale;
    m3scale = mplpot::m3scale;
    m4scale = mplpot::m4scale;
    m5scale = mplpot::m5scale;
  }
}

void empole(int vers) {
  if (empole_electyp == elec_t::coulomb)
    empole_coulomb(vers);
  else if (empole_electyp == elec_t::ewald)
    empole_ewald(vers);
}
TINKER_NAMESPACE_END
