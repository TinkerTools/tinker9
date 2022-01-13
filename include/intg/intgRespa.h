#pragma once
#include "baroBasic.h"
#include "enum.h"
#include "intgBasic.h"
#include "thermoBasic.h"

namespace tinker {
class RespaIntegrator : public BasicIntegrator
{
private:
   BasicThermostat* m_thermo;
   BasicBarostat* m_baro;
   int m_nrespa;

public:
   ~RespaIntegrator();
   RespaIntegrator(ThermostatEnum, BarostatEnum);
   void printDetail(FILE*) override;
   void kickoff() override;
   void dynamic(int, time_prec) override;

   static void updateRespaVelocity(time_prec tfast, time_prec tslow);
};
}
