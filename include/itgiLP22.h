#pragma once
#include "itgbBasic.h"
#include "itgiBasic.h"
#include "itgtNhc96.h"

namespace tinker {
class LP22Thermostat : public BasicThermostat
{
protected:
   Nhc96Thermostat* m_tpart;
   Nhc96Thermostat* m_tbaro;

public:
   static double* vbarKinetic();
   static void scaleVbarVelocity(double scale);

   LP22Thermostat();
   ~LP22Thermostat();
   void printDetail(FILE*) override;
   void control1b(time_prec timeStep, bool applyBaro);
   void control2b(time_prec timeStep, bool applyBaro);
};

class LP22Barostat : public BasicBarostat
{
protected:
   static constexpr int m_shapeArray[6][2] = {{0, 0}, {1, 1}, {2, 2},
                                             {0, 2}, {0, 1}, {1, 2}};
   double m_dofP;
   double m_fric;
   double m_rdn[3][3];
   int m_arrlen;

   void control12Impl(time_prec);

public:
   LP22Barostat();
   double dof() const;
   void printDetail(FILE*) override;
   void control1(time_prec timeStep) override;
   void control2(time_prec timeStep) override;
   void control3(time_prec timeStep) override;
};

// class LP22Integrator : public BasicIntegrator
// {};
}
