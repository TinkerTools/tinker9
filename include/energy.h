#pragma once
#include "time_scale.h"


#include "eangle.h"
#include "ebond.h"
#include "eopbend.h"
#include "epitors.h"
#include "estrbnd.h"
#include "etors.h"
#include "etortor.h"
#include "eurey.h"


#include "egeom.h"


#include "echarge.h"
#include "emplar.h"
#include "empole.h"
#include "epolar.h"
#include "evdw.h"


#include "echgtrn.h"
#include "ehippo.h"


namespace tinker {
/**
 * \ingroup mdegv
 * \brief Evaluate potential energy.
 * \param vers      Flag to select the version of energy routines.
 * \param tsflag    Time scale flag, a 32-bit encrypted integer.
 * \param tsconfig  Constant reference to a TimeScaleConfig object.
 *
 * \see TimeScaleConfig
 */
void energy_core(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);


/**
 * \ingroup mdegv
 * \brief First, energy buffers, virial buffers, gradient arrays, and count
 * buffers are set to 0. Then, evaluate energies, gradients, virials, and count
 * interactions. Last, update the global energy and virial tensor variables.
 * Counts are not updated until #count_reduce() is explicitly called. May skip
 * some steps based on the version parameter.
 */
void energy(int vers);


/**
 * \ingroup mdegv
 * \brief Zero-initialize and evaluate potential energy with time scale
 * configuration.
 * \see TimeScaleConfig
 */
void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);


void energy_data(rc_op);
}
