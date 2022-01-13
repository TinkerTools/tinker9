#pragma once
#include "baroBasic.h"
#include "enum.h"
#include "intgBasic.h"
#include "thermoBasic.h"

namespace tinker {
class VerletIntegrator : public BasicIntegrator
{
private:
   BasicThermostat* m_thermo;
   BasicBarostat* m_baro;

public:
   ~VerletIntegrator();
   VerletIntegrator(ThermostatEnum, BarostatEnum);
   void printDetail(FILE*) override;
   void kickoff() override;
   void dynamic(int, time_prec) override;
};
}
