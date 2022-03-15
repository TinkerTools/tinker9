#pragma once
#include "itgbLogVIso.h"
#include "itgtNhc06.h"

namespace tinker {
class Nhc06Barostat : public BasicBarostat
{
protected:
   Nhc06Thermostat* m_thermo;
   LogVIsoBarostat* m_baro;

public:
   Nhc06Barostat();
   ~Nhc06Barostat();
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control1(time_prec timeStep) override;
   void control2(time_prec timeStep) override;
   void control3(time_prec timeStep) override;
};
}
