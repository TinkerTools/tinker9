#pragma once
#include "glob.box.h"
#include "tool/rc_man.h"

namespace tinker {
void boxData(rc_op);
void boxExtent(double new_extent);
void boxSetDefault(const Box&);
void boxGetDefault(Box&);
void boxSetDefaultRecip();
void boxSetTinkerModule(const Box&);
void boxGetAxesAngles(const Box&, double& a, double& b, double& c, //
   double& alpha, double& beta, double& gamma);

/// \ingroup box
/// \brief Similar to Tinker `lattice` subroutine.
void boxLattice(Box& p, BoxShape sh, double a, double b, double c, //
   double alpha_deg, double beta_deg, double gamma_deg);

void boxCopyin();

/// \ingroup box
/// \brief Get the volume of the PBC box.
/// \note This function may calculate the volume on-the-fly instead of using
/// an internal variable to save the volume.
/// \note The volume is undefined for the unbound box.
real boxVolume();
}
