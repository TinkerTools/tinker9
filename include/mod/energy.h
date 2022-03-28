#pragma once
#include "ff/energybuffer.h"

namespace tinker {
TINKER_EXTERN int nelem_buffer;

/// \ingroup mdegv
/// \page mdegv
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

/// \ingroup mdegv
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

/// \var eksum
/// Kinetic energy.
/// \var ekin
/// Kinetic energy tensor.
TINKER_EXTERN energy_prec eksum;
TINKER_EXTERN energy_prec ekin[3][3];

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
