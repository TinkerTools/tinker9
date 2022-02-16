#include "itgbBerendsen.h"
#include "itgEnum.h"
#include "mdpt.h"

namespace tinker {
BarostatEnum BerendsenBarostat::getBarostatEnum() const
{
   return BarostatEnum::Berendsen;
}

void BerendsenBarostat::control2(time_prec dt)
{
   if (not applyBaro)
      return;

   berendsen_barostat(dt);
}
}
