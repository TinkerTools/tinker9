#include "md/integrator.h"
#include "tool/argkey.h"
#include "tool/ioprint.h"
#include <tinker/detail/bath.hh>

namespace tinker {
void BasicBarostat::printBasic(FILE* o)
{
   print(o, " Pressure           %12.1lf Atm\n", bath::atmsph);
   print(o, " Tau-Pressure       %12.1lf ps\n", bath::taupres);
   print(o, " Volume Trial       %12d\n", m_nbaro);
   if (semiiso)
      print(o, " Semiisotropic Fluctuation\n");
   else if (aniso)
      print(o, " Anisotropic Fluctuation\n");
   else
      print(o, " Isotropic Fluctuation\n");
}

BasicBarostat::BasicBarostat()
   : m_nbaro(1)
{
   int msave = m_nbaro;
   getKV("VOLUME-TRIAL", m_nbaro, msave);
}

BasicBarostat::~BasicBarostat() {}

void BasicBarostat::printDetail(FILE* o) {}

BarostatEnum BasicBarostat::getBarostatEnum() const
{
   return BarostatEnum::NONE;
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

BasicBarostat* BasicBarostat::create(BarostatEnum be)
{
   BasicBarostat* b = nullptr;
   switch (be) {
   case BarostatEnum::BERENDSEN:
      b = new BerendsenBarostat;
      break;
   case BarostatEnum::LP2022:
      b = new LP22Barostat;
      break;
   case BarostatEnum::MONTECARLO:
      b = new MonteCarloBarostat;
      break;
   case BarostatEnum::NHC2006:
      b = new Nhc06Barostat;
      break;
   default:
      b = new BasicBarostat;
      break;
   }
   return b;
}
}
