#pragma once
#include "itgbBasic.h"
#include "itgiBasic.h"
#include "itgpVerlet.h"
#include "itgtBasic.h"

namespace tinker {
enum class PropagatorEnum;
enum class ThermostatEnum;
enum class BarostatEnum;

class VerletIntegrator : public BasicIntegrator
{
protected:
   BasicThermostat* m_thermo;
   BasicBarostat* m_baro;
   int m_nrespa;
   void kickoff() override;

public:
   ~VerletIntegrator();
   VerletIntegrator(ThermostatEnum, BarostatEnum);
   VerletIntegrator(PropagatorEnum, ThermostatEnum, BarostatEnum) {}
   void printDetail(FILE*) override;
   void dynamic(int, time_prec) override;

   static void updateVelocity2(time_prec tfast, time_prec tslow);
};
}
