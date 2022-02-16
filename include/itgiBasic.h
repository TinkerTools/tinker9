#pragma once
#include "itgbBasic.h"
#include "itgpBasic.h"
#include "itgtBasic.h"

namespace tinker {
class BasicIntegrator : public IntegratorStaticData
{
protected:
   BasicPropagator* m_prop;
   BasicThermostat* m_thermo;
   BasicBarostat* m_baro;

   int vers1;
   bool save;

   void printBasic(FILE*);
   virtual void plan(int istep);

   virtual void kickoff() = 0;

public:
   BasicIntegrator(PropagatorEnum pe, ThermostatEnum te, BarostatEnum be);
   BasicIntegrator();
   virtual ~BasicIntegrator();
   virtual void printDetail(FILE*);
   virtual void dynamic(int istep, time_prec dt);
};
}
