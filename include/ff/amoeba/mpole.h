#pragma once

namespace tinker {
enum
{
   MPL_PME_0 = 0,
   MPL_PME_X = 1,
   MPL_PME_Y = 2,
   MPL_PME_Z = 3,
   MPL_PME_XX = 4,
   MPL_PME_YY = 5,
   MPL_PME_ZZ = 6,
   MPL_PME_XY = 7,
   MPL_PME_XZ = 8,
   MPL_PME_YZ = 9,
   MPL_TOTAL = 10,
   MPL_PME_YX = MPL_PME_XY,
   MPL_PME_ZX = MPL_PME_XZ,
   MPL_PME_ZY = MPL_PME_YZ,

   LFRM_NONE = 0,
   LFRM_Z_ONLY = 1,
   LFRM_Z_THEN_X = 2,
   LFRM_BISECTOR = 3,
   LFRM_Z_BISECT = 4,
   LFRM_3_FOLD = 5
};

/// \brief Local axis type and x,y,z-axis defining atoms for each multipole site.
struct LocalFrame
{
   int zaxis;  ///< Z-axis defining atom, starting from 0.
   int xaxis;  ///< X-axis defining atom, starting from 0.
   int yaxis;  ///< Y-axis defining atom, starting from ONE.
   int polaxe; ///< Local frame definition.
};

enum class UPred
{
   NONE,
   GEAR,
   ASPC,
   LSQR
};
}
