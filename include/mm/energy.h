#pragma once
#include "time_scale.h"

#include "amoeba/emplar.h"
#include "amoeba/empole.h"
#include "amoeba/epolar.h"
#include "evalence.h"
#include "evdw.h"
#include "pchg/echarge.h"
#include "pchg/echglj.h"

#include "hippo/echgtrn.h"
#include "hippo/edisp.h"
#include "hippo/ehippo.h"
#include "hippo/empole_chgpen.h"
#include "hippo/epolar_chgpen.h"
#include "hippo/erepel.h"

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
extern bool ecore_val;
extern bool ecore_vdw;
extern bool ecore_ele;

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

void energy_data(RcOp);
bool use_energi_vdw();
bool use_energi_elec();
}
