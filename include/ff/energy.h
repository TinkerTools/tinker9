#pragma once
#include "ff/atom.h"
#include "ff/box.h"
#include "mod/energy.h"
#include <map>
#include <string>

namespace tinker {
/// \ingroup mdegv
/// \brief Time scale configuration that assigns a group number to every energy
/// term. Up to 32 different group numbers are supported, from 0 to 31.
///
/// This class must:
///    - implement a C++ style `.at(arg)` member function where `arg` is the name
///      of an energy term;
///    - return the group to which `arg` belongs;
///    - throw an exception if `arg` is an invalid energy name.
///
/// If the k-th bit of the flag is set, all of the energies in `group k` will be
/// calculated.
///
/// #### Example 1.
/// Every energy term in the Verlet integrator is computed at every time-step.
/// Therefore, all of the terms are in the same group. This group (G) can be any
/// number from 0 to 31. Accordingly, the flag must be \f$ 2^G \f$.
///
/// #### Example 2.
/// There are usually two groups, "high frequency" (fast) and "low frequency"
/// (slow), of energy terms in the RESPA integrator. If one decides to assign 3
/// and 5 to these two groups, respectively, `tsconfig.at("ebond")` shall return
/// 3 and `tsconfig.at("evdw")` shall return 5.
///
/// The values of flags shall be:
///    - 8 (\f$ 2^3 \f$) for only the fast terms;
///    - 32 (\f$ 2^5 \f$) for only the slow terms;
///    - 40 (8 + 32) for both groups.
///
/// \note To use a new energy term in an existing integrator, the force field
/// developers are responsible to update the time scale configuration of this
/// integrator.
///
/// \todo Add or modify time scale configurations by keywords, e.g.
/// `"TIME-SCALE RESPA VDWTERM 0"`.
using TimeScaleConfig = std::map<std::string, int>;

/// \ingroup mdegv
/// \brief The default #TimeScaleConfig object. All of the energy terms are in
/// group 0, and can only be used via flag 1.
const TimeScaleConfig& defaultTSConfig();

/// \ingroup mdegv
/// \brief Evaluate potential energy.
/// \param vers      Flag to select the version of energy routines.
/// \param tsflag    Time scale flag, a 32-bit encrypted integer.
/// \param tsconfig  Constant reference to a TimeScaleConfig object.
///
/// \see TimeScaleConfig
void energy_core(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);

/// \ingroup mdegv
/// \brief First, energy buffers, virial buffers, gradient arrays, and count
/// buffers are set to 0. Then, evaluate energies, gradients, virials, and count
/// interactions. Last, update the global energy and virial tensor variables.
/// Counts are not updated until #countReduce() is explicitly called. May skip
/// some steps based on the version parameter.
void energy(int vers);

/// \ingroup mdegv
/// \brief Zero-initialize and evaluate potential energy with time scale
/// configuration.
/// \see TimeScaleConfig
void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);

void energyData(RcOp);
bool useEnergyVdw();
bool useEnergyElec();

/**
 * \ingroup mdegv
 * \brief
 * Zero out all of the counts, energies, gradients, and virials on device.
 */
void zero_egv(int vers);
/**
 * \ingroup mdegv
 * \brief
 * Zero out all of the counts, energies, gradients, and virials on device, with
 * parameter #rc_flag & calc::vmask.
 *
 * \see rc_flag
 * \see calc::vmask
 */
void zero_egv();

/**
 * \ingroup mdegv
 */
void scale_gradient(double scale, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z);
void sum_gradient(grad_prec* g0x, grad_prec* g0y, grad_prec* g0z, const grad_prec* g1x,
   const grad_prec* g1y, const grad_prec* g1z);
void sum_gradient(double scale, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z,
   const grad_prec* g1x, const grad_prec* g1y, const grad_prec* g1z);
void scale_gradient_acc(double, grad_prec*, grad_prec*, grad_prec*);
void sum_gradient_acc(
   grad_prec*, grad_prec*, grad_prec*, const grad_prec*, const grad_prec*, const grad_prec*);
void sum_gradient_acc(double, grad_prec*, grad_prec*, grad_prec*, const grad_prec*,
   const grad_prec*, const grad_prec*);

void copy_gradient(int vers, double* grdx, double* grdy, double* grdz, const grad_prec* gx_src,
   const grad_prec* gy_src, const grad_prec* gz_src, int queue);
/**
 * \ingroup mdegv
 * \brief
 * Copy total potential energy to another variable. Avoid accessing #esum
 * directly because it is not always available in the calculation, in which
 * case the output variable will not be changed.
 */
void copy_energy(int vers, energy_prec* eng);
/**
 * \ingroup mdegv
 * \brief
 * Copy the energy gradients from device to host.
 */
void copy_gradient(int vers, double* grdx, double* grdy, double* grdz, const grad_prec* gx_src,
   const grad_prec* gy_src, const grad_prec* gz_src);
/**
 * \ingroup mdegv
 * \brief
 * Copy the energy gradients from #gx, #gy, #gz to host.
 *
 * \see gx
 * \see gy
 * \see gz
 */
void copy_gradient(int vers, double* grdx, double* grdy, double* grdz);

void egvData(RcOp);
}
