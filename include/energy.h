#pragma once
#include "time_scale.h"


#include "e_angle.h"
#include "e_bond.h"
#include "e_opbend.h"
#include "e_pitors.h"
#include "e_strbnd.h"
#include "e_tors.h"
#include "e_tortor.h"
#include "e_urey.h"


#include "e_geom.h"


#include "e_mplar.h"
#include "e_mpole.h"
#include "e_polar.h"
#include "evdw.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup md_egv
 * \brief Evaluate the potential energy.
 *
 * First, energy buffers, virial buffers, #gx, #gy, #gz arrays, and count
 * buffers are set to 0.
 * Then, evaluate energies, gradients, virials, and count interactions.
 * Last, update the global variables #esum and #vir[9]; counts are not updated
 * until #count_reduce() is explicitly called.
 * May skip some steps based on the version parameter.
 *
 * \param vers      Flag to select the version of energy routines.
 * \param tsflag    Time scale flag, a 32-bit encrypted integer.
 * \param tsconfig  Constant reference to a TimeScaleConfig object.
 *
 * \see TimeScaleConfig
 */
void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);
void energy(int vers);


void energy_data(rc_op);
TINKER_NAMESPACE_END
