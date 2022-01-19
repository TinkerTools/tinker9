#include "energy.h"
#include "intg/baroMonteCarlo.h"
#include "intg/intgVerlet.h"
#include "mdcalc.h"
#include "mdpq.h"

namespace tinker {
VerletIntegrator::~VerletIntegrator()
{
   delete m_thermo;
   delete m_baro;
}

VerletIntegrator::VerletIntegrator(ThermostatEnum te, BarostatEnum be)
   : BasicIntegrator()
   , m_thermo(create(te))
   , m_baro(create(be))
{}

void VerletIntegrator::printDetail(FILE* o)
{
   m_thermo->printDetail(o);
   m_baro->printDetail(o);
   printBasic(o);
}

void VerletIntegrator::kickoff()
{
   energy(calc::v1);
}

void VerletIntegrator::dynamic(int istep, time_prec dt)
{
   int vers0 = rc_flag & calc::vmask;
   int vers1 = vers0;

   bool save = this->ifSave(istep);
   bool applyBaro = m_baro->ifApply(istep);
   bool mcbaro = false;
   if (applyBaro)
      if (dynamic_cast<MonteCarloBarostat*>(m_baro))
         mcbaro = true;
   // toggle off the calc::virial bit if Monte Carlo Barostat is in use
   if (mcbaro)
      vers1 &= ~calc::virial;
   // toggle off the calc::energy bit if neither save nor mcbaro
   if (not save and not mcbaro)
      vers1 &= ~calc::energy;

   time_prec dt2 = 0.5 * dt;

   if (applyBaro)
      m_baro->control1(dt);
   m_thermo->control1(dt);

   // v += a * dt/2
   updateVelocity(dt2);

   this->rattleSave();

   m_thermo->control3(dt);
   if (applyBaro)
      m_baro->control3(dt);

   updatePosition(dt);
   this->rattle(dt);
   copy_pos_to_xyz(true);

   energy(vers1);

   m_thermo->control4(dt);
   if (applyBaro)
      m_baro->control4(dt);

   // v += a * dt/2 and rattle2
   updateVelocity(dt2);
   this->rattle2(dt, vers1 & calc::virial);

   // full-step corrections
   m_thermo->control2(dt, save);
   if (applyBaro)
      m_baro->control2(dt);
}
}
