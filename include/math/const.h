#pragma once
#include "precision.h"
#include <cmath>

namespace tinker {
/// \ingroup math
/// \brief \f$ \pi \f$
constexpr real pi = M_PI;
/// \ingroup math
/// \brief \f$ \sqrt{\pi} \f$
constexpr real sqrtpi = 1.77245385090551602730;
/// \ingroup math
/// \brief \f$ 180/\pi \f$
constexpr real radian = 57.2957795130823208768;
/// \ingroup math
/// \brief \f$ \pi/180 \f$
constexpr real _1radian = 0.01745329251994329576924;
}
