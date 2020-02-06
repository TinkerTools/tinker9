#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
namespace calc {
/// \ingroup md
/// Use coordinates.
constexpr int xyz = 0x001;
/// \ingroup md
/// Use velocities.
constexpr int vel = 0x002;
/// \ingroup md
/// Use mass.
constexpr int mass = 0x004;
/// \ingroup md
/// Use multi-frame trajectory.
constexpr int traj = 0x008;


/// \ingroup md
/// Evaluate energy.
constexpr int energy = 0x010;
/// \ingroup md
/// Evaluate energy gradient.
constexpr int grad = 0x020;
/// \ingroup md
/// Evaluate virial tensor.
constexpr int virial = 0x040;
/// \ingroup md
/// Evaluate number of interactions.
constexpr int analyz = 0x080;


/// \ingroup md
/// Bits mask to clear energy-irrelevant flags.
constexpr int vmask = energy + grad + virial + analyz;
/// \ingroup md
/// Similar to basic Tinker energy routines.
/// Energy only.
constexpr int v0 = energy;
/// \ingroup md
/// Similar to version 1 Tinker energy routines.
/// Energy, gradient, and virial.
constexpr int v1 = energy + grad + virial;
/// \ingroup md
/// Similar to version 3 Tinker energy routines.
/// Energy and number of interactions.
constexpr int v3 = energy + analyz;
/// \ingroup md
/// Energy and gradient.
constexpr int v4 = energy + grad;
/// \ingroup md
/// Gradient only.
constexpr int v5 = grad;
/// \ingroup md
/// Gradient and virial.
constexpr int v6 = grad + virial;


/// \ingroup md
/// Run MD simulation.
constexpr int md = 0x100;
}


template <int USE>
struct EnergyVersion
{
   static constexpr int e = USE & calc::energy;
   static constexpr int a = USE & calc::analyz;
   static constexpr int g = USE & calc::grad;
   static constexpr int v = USE & calc::virial;
   static_assert(v ? g : true, "If calc::virial, must calc::grad.");
   static_assert(a ? e : true, "If calc::analyz, must calc::energy.");
};
TINKER_NAMESPACE_END
