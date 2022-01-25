#include "itgbBasic.h"
#include "itgEnum.h"
#include "tool/io_print.h"
#include <tinker/detail/bath.hh>

namespace tinker {
void BasicBarostat::printBasic(FILE* o)
{
   print(o, "\n");
   print(o, " Pressure         : %12.1lf Atm\n", bath::atmsph);
   print(o, " Tau-Pressure     : %12.1lf ps\n", bath::taupres);
   print(o, " NBaro            : %12d\n", m_nbaro);
}

BasicBarostat::BasicBarostat(BarostatEnum be)
   : m_baroEnum(be)
   , m_nbaro(1)
{}

BasicBarostat::BasicBarostat()
   : m_baroEnum(BarostatEnum::Null)
   , m_nbaro(1)
{}

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

BarostatEnum BasicBarostat::getBarostatEnum() const
{
   return m_baroEnum;
}
}
