#pragma once
#include "tool/energybuffer.h"

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

// mpole
namespace tinker {
TINKER_EXTERN LocalFrame* zaxis;
TINKER_EXTERN real (*pole)[MPL_TOTAL];
TINKER_EXTERN real (*rpole)[MPL_TOTAL];
TINKER_EXTERN real* trqx;
TINKER_EXTERN real* trqy;
TINKER_EXTERN real* trqz;

TINKER_EXTERN VirialBuffer vir_trq;

TINKER_EXTERN CountBuffer nem;
TINKER_EXTERN EnergyBuffer em;
TINKER_EXTERN VirialBuffer vir_em;
TINKER_EXTERN grad_prec* demx;
TINKER_EXTERN grad_prec* demy;
TINKER_EXTERN grad_prec* demz;
TINKER_EXTERN energy_prec energy_em;
TINKER_EXTERN virial_prec virial_em[9];
}

// mplpot
namespace tinker {
TINKER_EXTERN real m2scale;
TINKER_EXTERN real m3scale;
TINKER_EXTERN real m4scale;
TINKER_EXTERN real m5scale;

TINKER_EXTERN int nmexclude;
TINKER_EXTERN int (*mexclude)[2];
TINKER_EXTERN real* mexclude_scale;
}

// polar
namespace tinker {
TINKER_EXTERN real* polarity;
TINKER_EXTERN real* thole;
TINKER_EXTERN real* pdamp;

TINKER_EXTERN real (*udir)[3];
TINKER_EXTERN real (*udirp)[3];
TINKER_EXTERN real (*uind)[3];
TINKER_EXTERN real (*uinp)[3];

TINKER_EXTERN CountBuffer nep;
TINKER_EXTERN EnergyBuffer ep;
TINKER_EXTERN VirialBuffer vir_ep;
TINKER_EXTERN grad_prec* depx;
TINKER_EXTERN grad_prec* depy;
TINKER_EXTERN grad_prec* depz;
TINKER_EXTERN energy_prec energy_ep;
TINKER_EXTERN virial_prec virial_ep[9];

TINKER_EXTERN real* polarity_inv;

TINKER_EXTERN real (*ufld)[3];
TINKER_EXTERN real (*dufld)[6];

TINKER_EXTERN real (*work01_)[3];
TINKER_EXTERN real (*work02_)[3];
TINKER_EXTERN real (*work03_)[3];
TINKER_EXTERN real (*work04_)[3];
TINKER_EXTERN real (*work05_)[3];
TINKER_EXTERN real (*work06_)[3];
TINKER_EXTERN real (*work07_)[3];
TINKER_EXTERN real (*work08_)[3];
TINKER_EXTERN real (*work09_)[3];
TINKER_EXTERN real (*work10_)[3];
}

// polpot
namespace tinker {
TINKER_EXTERN real u1scale;
TINKER_EXTERN real u2scale;
TINKER_EXTERN real u3scale;
TINKER_EXTERN real u4scale;
TINKER_EXTERN real d1scale;
TINKER_EXTERN real d2scale;
TINKER_EXTERN real d3scale;
TINKER_EXTERN real d4scale;
TINKER_EXTERN real p2scale;
TINKER_EXTERN real p3scale;
TINKER_EXTERN real p4scale;
TINKER_EXTERN real p5scale;
TINKER_EXTERN real p2iscale;
TINKER_EXTERN real p3iscale;
TINKER_EXTERN real p4iscale;
TINKER_EXTERN real p5iscale;

TINKER_EXTERN real udiag;

TINKER_EXTERN int nuexclude;
TINKER_EXTERN int (*uexclude)[2];
TINKER_EXTERN real* uexclude_scale;

TINKER_EXTERN int ndpexclude;
TINKER_EXTERN int (*dpexclude)[2];
TINKER_EXTERN real (*dpexclude_scale)[2];

TINKER_EXTERN int ndpuexclude;
TINKER_EXTERN int (*dpuexclude)[2];
TINKER_EXTERN real (*dpuexclude_scale)[3];
}

// mplar
namespace tinker {
TINKER_EXTERN int nmdpuexclude;
TINKER_EXTERN int (*mdpuexclude)[2];
TINKER_EXTERN real (*mdpuexclude_scale)[4];
}

// uprior
namespace tinker {
TINKER_EXTERN UPred polpred;
TINKER_EXTERN int maxualt;
TINKER_EXTERN int nualt;
TINKER_EXTERN real (*udalt_00)[3];
TINKER_EXTERN real (*udalt_01)[3];
TINKER_EXTERN real (*udalt_02)[3];
TINKER_EXTERN real (*udalt_03)[3];
TINKER_EXTERN real (*udalt_04)[3];
TINKER_EXTERN real (*udalt_05)[3];
TINKER_EXTERN real (*udalt_06)[3];
TINKER_EXTERN real (*udalt_07)[3];
TINKER_EXTERN real (*udalt_08)[3];
TINKER_EXTERN real (*udalt_09)[3];
TINKER_EXTERN real (*udalt_10)[3];
TINKER_EXTERN real (*udalt_11)[3];
TINKER_EXTERN real (*udalt_12)[3];
TINKER_EXTERN real (*udalt_13)[3];
TINKER_EXTERN real (*udalt_14)[3];
TINKER_EXTERN real (*udalt_15)[3];
TINKER_EXTERN real (*upalt_00)[3];
TINKER_EXTERN real (*upalt_01)[3];
TINKER_EXTERN real (*upalt_02)[3];
TINKER_EXTERN real (*upalt_03)[3];
TINKER_EXTERN real (*upalt_04)[3];
TINKER_EXTERN real (*upalt_05)[3];
TINKER_EXTERN real (*upalt_06)[3];
TINKER_EXTERN real (*upalt_07)[3];
TINKER_EXTERN real (*upalt_08)[3];
TINKER_EXTERN real (*upalt_09)[3];
TINKER_EXTERN real (*upalt_10)[3];
TINKER_EXTERN real (*upalt_11)[3];
TINKER_EXTERN real (*upalt_12)[3];
TINKER_EXTERN real (*upalt_13)[3];
TINKER_EXTERN real (*upalt_14)[3];
TINKER_EXTERN real (*upalt_15)[3];

TINKER_EXTERN real* udalt_lsqr_a;
TINKER_EXTERN real* upalt_lsqr_a;
TINKER_EXTERN real* udalt_lsqr_b;
TINKER_EXTERN real* upalt_lsqr_b;
}
