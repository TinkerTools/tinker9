#pragma once
#include "itgbBasic.h"
#include "itgiBasic.h"
#include "itgtNhc.h"

namespace tinker {
class Nhc06Thermostat : public BasicThermostat
{
protected:
   NhcThermostat* m_tpart;
   NhcThermostat* m_tbaro;

public:
   static double* vbarKinetic();
   static void scaleVbarVelocity(double scale);

   Nhc06Thermostat();
   ~Nhc06Thermostat();
   void printDetail(FILE*) override;
   void control1b(time_prec timeStep, bool applyBaro);
   void control2b(time_prec timeStep, bool applyBaro);
};

class Nhc06Barostat : public BasicBarostat
{
protected:
   double m_dofP;
   void control12Impl(time_prec timeStep);

public:
   Nhc06Barostat();
   double dof() const;
   void printDetail(FILE*) override;
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
   Nhc06Integrator(bool useVerlet) {}
   Nhc06Integrator();
   ~Nhc06Integrator();
   void printDetail(FILE*) override;
   void dynamic(int, time_prec) override;

   static void updatePositionPedantic(time_prec t);
   static void updateVelocityPedantic(time_prec t, double velbar);
};
}
