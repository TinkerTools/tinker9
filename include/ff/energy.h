#pragma once
#include "ff/atom.h"
#include "ff/box.h"
#include "ff/egvop.h"
#include "ff/energybuffer.h"
#include "ff/timescale.h"
#include "tool/darray.h"

namespace tinker {
/// \ingroup egv
/// \brief Evaluate potential energy.
/// \param vers      Flag to select the version of energy routines.
/// \param tsflag    Time scale flag, a 32-bit encrypted integer.
/// \param tsconfig  Constant reference to a TimeScaleConfig object.
///
/// \see TimeScaleConfig
void energy_core(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);

/// \ingroup egv
/// \brief First, energy buffers, virial buffers, gradient arrays, and count
/// buffers are set to 0. Then, evaluate energies, gradients, virials, and count
/// interactions. Last, update the global energy and virial tensor variables.
/// Counts are not updated until #countReduce() is explicitly called. May skip
/// some steps based on the version parameter.
void energy(int vers);

/// \ingroup egv
/// \brief Zero-initialize and evaluate potential energy with time scale configuration.
/// \see TimeScaleConfig
void energy(int vers, unsigned tsflag, const TimeScaleConfig& tsconfig);

/// \ingroup egv
void energyData(RcOp);

/// \ingroup egv
bool useEnergyVdw();

/// \ingroup egv
bool useEnergyElec();

/// \ingroup egv
void egvData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup egv
/// \page egv
///
/// Charge, Multipole, Polarization, and HIPPO Charge Transfer Terms
///
/// |           | .not. analyz     | analyz    |
/// |-----------|------------------|-----------|
/// | ebuf_elec | new array        | nullptr   |
/// | vbuf_elec | new array        | nullptr   |
/// | g_elec    | new array        | new array |
/// | em        | alias: ebuf_elec | new array |
/// | vir_em    | alias: vbuf_elec | new array |
/// | dem       | alias: g_elec    | new array |
///
/// VDW, HIPPO Dispersion, and HIPPO Repulsion Terms
///
/// |          | .not. analyz    | analyz    |
/// |----------|-----------------|-----------|
/// | ebuf_vdw | new array       | nullptr   |
/// | vbuf_vdw | new array       | nullptr   |
/// | g_vdw    | new array       | new array |
/// | ev       | alias: ebuf_vdw | new array |
/// | vir_ev   | alias: vbuf_vdw | new array |
/// | dev      | alias: g_vdw    | new array |
///
/// Valence Terms
///
/// |        | .not. analyz | analyz    |
/// |--------|--------------|-----------|
/// | ebuf   | new array    | nullptr   |
/// | vbuf   | new array    | nullptr   |
/// | g      | new array    | new array |
/// | eb     | alias: ebuf  | new array |
/// | vir_eb | alias: vbuf  | new array |
/// | deb    | alias: g     | new array |

/// \ingroup egv
/// \{
/// \var gx
/// \brief Gradient of valence terms; also works as total gradient array.
/// \var gy
/// \copydoc gx
/// \var gz
/// \copydoc gx
TINKER_EXTERN grad_prec* gx;
TINKER_EXTERN grad_prec* gy;
TINKER_EXTERN grad_prec* gz;

/// \var gx_vdw
/// Gradient of VDW (HIPPO dispersion, HIPPO repulsion) terms.
/// \var gy_vdw
/// \copydoc gx_vdw
/// \var gz_vdw
/// \copydoc gx_vdw
TINKER_EXTERN grad_prec* gx_vdw;
TINKER_EXTERN grad_prec* gy_vdw;
TINKER_EXTERN grad_prec* gz_vdw;

/// \var gx_elec
/// Gradient of electrostatic (charge, multipole, polarization, HIPPO charge transfer) terms.
/// \var gy_elec
/// \copydoc gx_elec
/// \var gz_elec
/// \copydoc gx_elec
TINKER_EXTERN grad_prec* gx_elec;
TINKER_EXTERN grad_prec* gy_elec;
TINKER_EXTERN grad_prec* gz_elec;

/// \var eng_buf
/// Energy buffer for the valence terms.
/// \var eng_buf_vdw
/// Energy buffer for the vdw terms.
/// \var eng_buf_elec
/// Energy buffer for the electrostatic terms.
TINKER_EXTERN EnergyBuffer eng_buf;
TINKER_EXTERN EnergyBuffer eng_buf_vdw;
TINKER_EXTERN EnergyBuffer eng_buf_elec;

/// \var vir_buf
/// Virial buffer for the valence terms.
/// \var vir_buf_vdw
/// Virial buffer for the vdw terms.
/// \var vir_buf_elec
/// Virial buffer for the electrostatic terms.
TINKER_EXTERN VirialBuffer vir_buf;
TINKER_EXTERN VirialBuffer vir_buf_vdw;
TINKER_EXTERN VirialBuffer vir_buf_elec;

/// \var esum
/// Total potential energy.
/// \var energy_valence
/// Total valence energy.
/// \var energy_vdw
/// Total vdw energy.
/// \var energy_elec
/// Total electrostatic energy.
TINKER_EXTERN energy_prec esum;
TINKER_EXTERN energy_prec energy_valence;
TINKER_EXTERN energy_prec energy_vdw;
TINKER_EXTERN energy_prec energy_elec;

/// \var vir
/// Total potential virial tensor.
/// \var virial_valence
/// Total valence virial tensor.
/// \var virial_vdw
/// Total vdw virial tensor.
/// \var virial_elec
/// Total electrostatic virial tensor.
TINKER_EXTERN virial_prec vir[9];
TINKER_EXTERN virial_prec virial_valence[9];
TINKER_EXTERN virial_prec virial_vdw[9];
TINKER_EXTERN virial_prec virial_elec[9];
/// \}
}
