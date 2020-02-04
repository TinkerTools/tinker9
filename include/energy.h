#pragma once


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
#include "e_vdw.h"


#include <map>


TINKER_NAMESPACE_BEGIN
using TimeScaleConfig = std::map<const char*, int>;
const TimeScaleConfig& default_tsconfig();


/**
 * \param time_scale  32-bit encrypted integer flag.
 * The integrator will assign a "group number" ranging from 0 to 31 to every
 * energy term. If the i-th bit of `time_scale` is set, all of the energy terms
 * of the i-th group will be calculated.
 *
 * For instance, 0x04 (0b0100) will only calculate group 2; 0x05 (0b0101) will
 * calculate both group 0 and group 2.
 *
 * \param tsconfig    Constant reference to a TimeScaleConfig object.
 * This object must implement a C++ style `.at(arg)` member function and throw
 * an exception if `arg` is an invalid name of the energy term.
 *
 * For instance, every term in the Verlet integrator is computed at every
 * time-step, therefore, all of the terms have the same group value (G). If G is
 * 0, the 0-th bit `time_scale` must be set; if G is 3, the 3-th bit
 * of `time_scale` must be set; etc.
 *
 * Another example is the RESPA integrator where we usually have two groups for
 * "high frequency" (fast) and "low frequency" (slow) terms. If one decides to
 * assign 3 and 4 to these two groups, respectively, `tsconfig.at("evdw")` must
 * return 4 and `tsconfig.at("ebond")` must return 3; to calculate the slow
 * terms, `time_scale` should be 0x08; to calculate the fast terms, `time_scale`
 * should be 0x04; to calculate both groups, `time_scale` should be 0x0C.
 *
 * \note To use a new energy term with an old integrator, the force field
 * developers are responsible to update the time scale configuration of this
 * integrator.
 */
void energy_potential(int vers, int time_scale = 1,
                      const TimeScaleConfig& tsconfig = default_tsconfig());
TINKER_NAMESPACE_END
