#pragma once
#include "itgbBasic.h"
#include "itgtNhc06.h"

namespace tinker {
class LP22Barostat : public BasicBarostat
{
protected:
   Nhc06Thermostat* m_thermo;
   BasicBarostat* m_baro;

public:
   LP22Barostat();
   ~LP22Barostat();
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control1(time_prec timeStep) override;
   void control2(time_prec timeStep) override;
   void control3(time_prec timeStep) override;
};
}
