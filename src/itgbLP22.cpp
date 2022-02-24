#include "itgbLP22.h"
#include "itgbLogVAniso.h"
#include "itgbLogVIso.h"
#include <tinker/detail/stodyn.hh>

namespace tinker {
LP22Barostat::LP22Barostat()
   : BasicBarostat()
   , m_baro(nullptr)
{
   if (aniso)
      m_baro = new LogVAnisoBarostat(stodyn::friction);
   else
      m_baro = new LogVIsoBarostat(stodyn::friction);
}

LP22Barostat::~LP22Barostat()
{
   delete m_baro;
}

void LP22Barostat::printDetail(FILE* o)
{
   m_baro->printDetail(o);
}

BarostatEnum LP22Barostat::getBarostatEnum() const
{
   return BarostatEnum::LP2022;
}

void LP22Barostat::control1(time_prec dt)
{
   __PlaceHolderMessage("LP22Barostat::control1T");
   m_baro->control1(dt);
}

void LP22Barostat::control2(time_prec dt)
{
   m_baro->control2(dt);
   __PlaceHolderMessage("LP22Barostat::control2T");
}

void LP22Barostat::control3(time_prec dt)
{
   m_baro->control3(dt);
}
}
