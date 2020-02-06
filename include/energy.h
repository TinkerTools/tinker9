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
#include "e_vdw.h"


TINKER_NAMESPACE_BEGIN
/**
 * \brief Evaluate the energy potential.
 *
 * First, energy buffers, virial buffers, and gx, gy, gz arrays are set to 0.
 * Then, evaluate energies, gradients, and virials.
 * Last, update the global variables esum and vir[9].
 * May skip some calculations based on the version parameter.
 *
 * \param vers      Flag to select the version of energy routines.
 *
 * \param tsflag    Time scale flag, a 32-bit encrypted integer.
 * The integrator will assign a "group number" ranging from 0 to 31 to every
 * energy term. If the i-th bit of `tsflag` is set, all of the energy terms
 * in the i-th group will be calculated.
 *
 * For instance, 4 (0b0100, 2**2) will only calculate group 2; 0x05 (0b0101,
 * 2**2 + 2**0) will calculate both group 2 and group 0.
 *
 * \param tsconfig  Constant reference to a TimeScaleConfig object.
 * This object must implement a C++ style `.at(arg)` member function and throw
 * an exception if `arg` is an invalid name of the energy term.
 *
 * For instance, every term in the Verlet integrator is computed at every
 * time-step, therefore, all of the terms have the same group value (G). If G is
 * 0, the 0-th bit `tsflag` must be set; if G is 3, the 3-th bit of `tsflag`
 * must be set; etc.
 *
 * Another example is the RESPA integrator where there are often two groups for
 * "high frequency" (fast) and "low frequency" (slow) terms. If one decides to
 * assign 3 and 5 to these two groups, respectively, `tsconfig.at("ebond")` must
 * return 3 and `tsconfig.at("evdw")` must return 5; to calculate the fast
 * terms, `tsflag` should be 8 (2**3); to calculate the slow terms, `tsflag`
 * should be 32 (2**5); to calculate both groups, `tsflag` should be 40 (2**3 +
 * 2**5).
 *
 * \note To use a new energy term in an old integrator, the force field
 * developers are responsible to update the time scale configuration of this
 * integrator.
 */
void energy_potential(int vers, int tsflag = 1,
                      const TimeScaleConfig& tsconfig = default_tsconfig());
/**
 * \brief If used and output pointer are not null: copy esum to eng;
 * copy gx, gy, gz to grdx, grdy, grdz, respectively;
 * copy vir[9] to virial.
 */
void copy_energy(int vers, real* eng, real* grdx, real* grdy, real* grdz,
                 real* virial);
TINKER_NAMESPACE_END
