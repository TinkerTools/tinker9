#include "md/integrator.h"
#include "md/misc.h"
#include "tool/ioprint.h"

namespace tinker {
void BerendsenBarostat::printDetail(FILE* o)
{
   print(o, "\n");
   print(o, " Berendsen Barostat\n");
   printBasic(o);
}

BarostatEnum BerendsenBarostat::getBarostatEnum() const
{
   return BarostatEnum::Berendsen;
}

void BerendsenBarostat::control2(time_prec dt)
{
   if (not applyBaro)
      return;

   berendsenBarostat(dt);
}
}
