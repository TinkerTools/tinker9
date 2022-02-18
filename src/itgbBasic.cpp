#include "itgbBasic.h"
#include "itgbMonteCarlo.h"
#include "tinker_rt.h"
#include "tool/io_print.h"
#include <tinker/detail/bath.hh>

namespace tinker {
void BasicBarostat::printBasic(FILE* o)
{
   print(o, " Pressure           %12.1lf Atm\n", bath::atmsph);
   print(o, " Tau-Pressure       %12.1lf ps\n", bath::taupres);
   print(o, " Volume Trial       %12d\n", m_nbaro);
}

BasicBarostat::BasicBarostat()
   : m_nbaro(1)
{
   int msave = m_nbaro;
   get_kv("VOLUME-TRIAL", m_nbaro, msave);
}

BasicBarostat::~BasicBarostat() {}

void BasicBarostat::printDetail(FILE* o) {}

BarostatEnum BasicBarostat::getBarostatEnum() const
{
   return BarostatEnum::Null;
}

bool BasicBarostat::ifApply(int istep)
{
   applyBaro = ((istep % m_nbaro) == 0);
   return applyBaro;
}

bool BasicBarostat::ifApply() const
{
   return applyBaro;
}
}
