#include "intg/baroBasic.h"
#include "tinker_rt.h"
#include "tool/io_print.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>

namespace tinker {
void BasicBarostat::printBasic(FILE* o)
{
   if (not inform::verbose)
      return;

   print(o, " Pressure         : %12.1lf\n", bath::atmsph);
   print(o, " Tau-Pressure     : %12.1lf\n", bath::taupres);
   print(o, " NBaro            : %12d\n", m_nbaro);
   print(o, "\n");
}

BasicBarostat::BasicBarostat()
{
   get_kv("VOLUME-TRIAL", m_nbaro, 1);
}

BasicBarostat::~BasicBarostat() {}

void BasicBarostat::printDetail(FILE* o)
{
   printBasic(o);
}

bool BasicBarostat::ifApply(int istep)
{
   int i = istep % m_nbaro;
   return i == 0;
}
}
