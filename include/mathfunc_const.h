#pragma once
#include "macro.h"
#include <cmath>


TINKER_NAMESPACE_BEGIN
/// \ingroup math
/// \f$ \sqrt[6]{2} \f$
constexpr real twosix = 1.12246204830937298143;
/// \ingroup math
/// \f$ \sqrt{2} \f$
constexpr real sqrttwo = 1.41421356237309504880;
/// \ingroup math
/// \f$ \sqrt{3} \f$
constexpr real sqrtthree = 1.73205080756887729353;


/// \ingroup math
/// \f$ exp(1) \f$
constexpr real elog = M_E;
/// \ingroup math
/// \f$ ln(10) \f$
constexpr real logten = M_LN10;


/// \ingroup math
/// \f$ \pi \f$
constexpr real pi = M_PI;
/// \ingroup math
/// \f$ \sqrt{\pi} \f$
constexpr real sqrtpi = 1.77245385090551602730;
/// \ingroup math
/// \f$ 180/\pi \f$
constexpr real radian = 57.2957795130823208768;
/// \ingroup math
/// \f$ \pi/180 \f$
constexpr real _1radian = 0.01745329251994329576924;
TINKER_NAMESPACE_END
