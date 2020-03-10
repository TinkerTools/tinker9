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


#include "emplar.h"
#include "empole.h"
#include "epolar.h"
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
