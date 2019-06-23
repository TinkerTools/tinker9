#include "gpu/decl_potent.h"
#include <cassert>
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
int use_potent(potent_t term) {
  int val = 0;
  switch (term) {
  case bond_term:
    val = potent::use_bond;
    break;
  case angle_term:
    val = potent::use_angle;
    break;
  default:
    assert(false);
    break;
  }
  return val;
}

int count_bonded_term(potent_t term) {
  int val = -1;
  if (use_potent(term)) {
    switch (term) {
    case bond_term:
      val = bndstr::nbond;
      break;
    case angle_term:
      val = angbnd::nangle;
      break;
    default:
      assert(false);
      break;
    }
  }
  return val;
}
}
TINKER_NAMESPACE_END
