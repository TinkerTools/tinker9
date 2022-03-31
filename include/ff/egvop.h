#pragma once
#include "precision.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup egv
/// \brief Zero out all of the counts, energies, gradients, and virials on device.
void zeroEGV(int vers);

/// \ingroup egv
/// \brief Zero out all of the counts, energies, gradients, and virials on device,
/// with extra information in #rc_flag.
///
/// \see rc_flag
/// \see calc::vmask
void zeroEGV();

/// \ingroup egv
/// \{
void scaleGradient(double scale, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z);

void sumGradient(grad_prec* g0x, grad_prec* g0y, grad_prec* g0z, //
   const grad_prec* g1x, const grad_prec* g1y, const grad_prec* g1z);

void sumGradient(double scale, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z, //
   const grad_prec* g1x, const grad_prec* g1y, const grad_prec* g1z);

void copyGradient(int vers, double* grdx, double* grdy, double* grdz, //
   const grad_prec* gx_src, const grad_prec* gy_src, const grad_prec* gz_src, int queue);

///  \brief Copy the energy gradients from device to host.
void copyGradient(int vers, double* grdx, double* grdy, double* grdz, //
   const grad_prec* gx_src, const grad_prec* gy_src, const grad_prec* gz_src);

/// \brief Copy the energy gradients from #gx, #gy, #gz to host.
///
/// \see gx
/// \see gy
/// \see gz
void copyGradient(int vers, double* grdx, double* grdy, double* grdz);

/// \brief Copy total potential energy to another variable. Avoid accessing #esum
/// directly because it is not always available in the calculation, in which
/// case the output variable will not be changed.
void copyEnergy(int vers, energy_prec* eng);
/// \}
}
