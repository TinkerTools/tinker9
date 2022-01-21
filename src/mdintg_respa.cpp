#include "mathfunc_pow2.h"
#include "mdintg.h"


namespace tinker {
grad_prec *gx1, *gy1, *gz1;
grad_prec *gx2, *gy2, *gz2;


const TimeScaleConfig& respa_tsconfig()
{
   constexpr int fast = floor_log2_constexpr(RESPA_FAST); // short-range
   constexpr int slow = floor_log2_constexpr(RESPA_SLOW); // long-range
   static TimeScaleConfig tsconfig{
      {"ebond", fast},         {"eangle", fast},        {"estrbnd", fast},
      {"eurey", fast},         {"eopbend", fast},       {"etors", fast},
      {"eimprop", fast},       {"eimptor", fast},       {"epitors", fast},
      {"estrtor", fast},       {"eangtor", fast},       {"etortor", fast},
      {"egeom", fast},

      {"evalence", fast},

      {"evdw", slow},

      {"echarge", slow},       {"echglj", slow},

      {"emplar", slow},        {"empole", slow},        {"epolar", slow},

      {"empole_chgpen", slow}, {"epolar_chgpen", slow},

      {"echgtrn", slow},       {"edisp", slow},         {"erepel", slow},
      {"ehippo", slow},
   };
   return tsconfig;
}
}
