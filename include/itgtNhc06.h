#pragma once
#include "itgtNhc.h"

namespace tinker {
class Nhc06Thermostat : public BasicThermostat
{
protected:
   NhcDevice* m_tpart;
   NhcDevice* m_tbaro;

public:
   Nhc06Thermostat();
   ~Nhc06Thermostat();

   void printDetail(FILE*) override;
   void control1(time_prec dt) override;
   void control2(time_prec dt, bool save) override;

   static double kineticRattleGroup();
   static void scaleVelocityRattleGroup(double scale);

   static double kineticVbar();
   static void scaleVelocityVbar(double scale);

   static double dofVbar();
};
}
