#ifndef TINKER_GPU_E_POTENTIAL_H_
#define TINKER_GPU_E_POTENTIAL_H_

#include "e_angle.h"
#include "e_bond.h"
#include "e_opbend.h"
#include "e_pitors.h"
#include "e_strbnd.h"
#include "e_tors.h"
#include "e_tortor.h"
#include "e_urey.h"

#include "e_mpole.h"
#include "e_polar.h"
#include "e_vdw.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(rc_t rc);

void energy_potential(int vers);
}
TINKER_NAMESPACE_END

#endif
