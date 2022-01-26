#pragma once
#include "itgbBasic.h"
#include "itgiBasic.h"
#include "itgtNhc96.h"

namespace tinker {
class Nhc06Thermostat : public BasicThermostat
{
protected:
   Nhc96Thermostat* m_tpart;
   Nhc96Thermostat* m_tbaro;

   static double* vbarKinetic();
   static void scaleVbarVelocity(double scale);

public:
   Nhc06Thermostat();
   ~Nhc06Thermostat();
   void control1(time_prec timeStep, bool applyBaro);
   void control2(time_prec timeStep, bool save, bool applyBaro);
   double dof() const;
};

class Nhc06Barostat : public BasicBarostat
{
protected:
   double m_dofP;
   void control12Impl(time_prec timeStep);

public:
   double dof() const;
   Nhc06Barostat();
   void control1(time_prec timeStep) override;
   void control2(time_prec timeStep) override;
   void control3(time_prec timeStep) override;
};

class Nhc06Integrator : public BasicIntegrator
{
protected:
   Nhc06Thermostat* m_thermo;
   Nhc06Barostat* m_baro;
   bool m_pedantic;
   void kickoff() override;

public:
   Nhc06Integrator();
   ~Nhc06Integrator();
   void printDetail(FILE*) override;
   void dynamic(int, time_prec) override;

   static void updatePositionPedantic(time_prec t);
   static void updateVelocityPedantic(time_prec t, double velbar);
};
}
