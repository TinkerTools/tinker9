#pragma once
#include <map>
#include <string>

namespace tinker {
/// \ingroup egv
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

/// \ingroup egv
/// \brief The default #TimeScaleConfig object. All of the energy terms are in
/// group 0, and can only be used via flag 1.
const TimeScaleConfig& defaultTSConfig();
}
