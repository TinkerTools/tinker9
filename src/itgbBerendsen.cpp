#include "itgbBerendsen.h"
#include "itgEnum.h"
#include "mdpt.h"
#include "tool/io_print.h"

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

   berendsen_barostat(dt);
}
}
