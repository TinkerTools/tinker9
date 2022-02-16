#include "itgbNhc06.h"
#include "itgbLogVAniso.h"
#include "itgbLogVIso.h"
#include <tinker/detail/bath.hh>

namespace tinker {
Nhc06Barostat::Nhc06Barostat()
   : BasicBarostat()
   , m_baro(nullptr)
{
   if (bath::anisotrop)
      m_baro = new LogVAnisoBarostat(0.0);
   else
      m_baro = new LogVIsoBarostat(0.0);
}

Nhc06Barostat::~Nhc06Barostat()
{
   delete m_baro;
}

void Nhc06Barostat::printDetail(FILE* o)
{
   m_baro->printDetail(o);
}

BarostatEnum Nhc06Barostat::getBarostatEnum() const
{
   return BarostatEnum::Nhc2006;
}

void Nhc06Barostat::control1(time_prec dt)
{
   __PlaceHolderMessage("Nhc06Barostat::control1T");
   m_baro->control1(dt);
}

void Nhc06Barostat::control2(time_prec dt)
{
   m_baro->control2(dt);
   __PlaceHolderMessage("Nhc06Barostat::control2T");
}

void Nhc06Barostat::control3(time_prec dt)
{
   m_baro->control3(dt);
}
}
