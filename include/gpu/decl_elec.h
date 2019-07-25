#ifndef TINKER_GPU_DECL_ELEC_H_
#define TINKER_GPU_DECL_ELEC_H_

#include "mod_elec.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/// @return 0 if no multipole electrostatics is involved; otherise, non-zero
int use_elec();

/**
 * @param op  construct or destruct multipole electrostatics
 *
 * has no effect if no multipole electrostatics is involved
 */
void elec_data(rc_t rc);

/**
 * initializes the electrostatics calculation
 *
 * Input:
 * @param vers  selects the code path
 *
 * Output:
 * zero torque (if used);
 * zero torque-related virial (if used);
 * call chkpole() and rotpole();
 * if use pme, initialize some pme data structures.
 */
void elec_init(int vers);

/**
 * takes the torque values on a single site defined by a local coordinate frame
 * and converts to Cartesian forces on the original site and sites specifying
 * the local frame.
 *
 * Input:
 * x, y, and z coordinates;
 * x, y, and z torques;
 * local frame definitions;
 * @param vers  selects the code path
 *
 * Output:
 * According to the value of @param ver, it may add add energy gradients, add
 * torque-related virial, or do nothing.
 */
void torque(int vers);
}
TINKER_NAMESPACE_END

#endif
