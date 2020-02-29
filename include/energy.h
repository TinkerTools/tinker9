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
 * \ingroup egv
 * \brief Evaluate the potential energy.
 *
 * First, energy buffers, virial buffers, #gx, #gy, #gz arrays, and count
 * buffers are set to 0.
 * Then, evaluate energies, gradients, virials, and count interactions.
 * Last, update the global variables #esum and #vir[9]; counts are not updated
 * until #get_count() is explicitly called.
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

/**
 * \ingroup egv
 * \brief Copy energy, gradients, and/or virials to other variables.
 *
 * If the output pointers are not null and are different from the sources:
 *    - copy #esum to `eng`;
 *    - non-blockingly copy #gx, #gy, #gz to `grdx`, `grdy`, and `grdz`,
 *      respectively;
 *    - copy #vir[9] to `virial[9]`.
 *
 * \param vers  Energy version flag to select data to be copied.
 */
void copy_energy(int vers, energy_prec* eng, real* grdx, real* grdy, real* grdz,
                 real* virial);
TINKER_NAMESPACE_END
