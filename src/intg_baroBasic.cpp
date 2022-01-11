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

int BasicBarostat::nbaro() const
{
   return m_nbaro;
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

void BasicBarostat::control(double) {}

void BasicBarostat::controlAfter(double) {}
}
