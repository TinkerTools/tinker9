#include "itgbNhc06.h"
#include "itgbLogVAniso.h"
#include "itgbLogVIso.h"
#include "tool/io.h"

namespace tinker {
Nhc06Barostat::Nhc06Barostat()
   : BasicBarostat()
   , m_thermo(new Nhc06Thermostat)
   , m_baro(nullptr)
{
   m_baro = new LogVIsoBarostat(0.0);
}

Nhc06Barostat::~Nhc06Barostat()
{
   delete m_baro;
   delete m_thermo;
}

void Nhc06Barostat::printDetail(FILE* o)
{
   print(o, "\n");
   m_baro->printDetail(o);
}

BarostatEnum Nhc06Barostat::getBarostatEnum() const
{
   return BarostatEnum::Nhc2006;
}

void Nhc06Barostat::control1(time_prec dt)
{
   m_thermo->control1(dt);
   m_baro->control1(dt);
}

void Nhc06Barostat::control2(time_prec dt)
{
   m_baro->control2(dt);
   m_thermo->control2(dt, false);
}

void Nhc06Barostat::control3(time_prec dt)
{
   m_baro->control3(dt);
}
}
