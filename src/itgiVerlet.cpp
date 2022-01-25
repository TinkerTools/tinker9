#include "itgiVerlet.h"
#include "energy.h"
#include "itgEnum.h"
#include "itgbMonteCarlo.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"

namespace tinker {
void VerletIntegrator::kickoff()
{
   energy(calc::grad);
}

VerletIntegrator::~VerletIntegrator()
{
   delete m_thermo;
   delete m_baro;
}

VerletIntegrator::VerletIntegrator(ThermostatEnum te, BarostatEnum be)
   : BasicIntegrator()
   , m_thermo(create(te))
   , m_baro(create(be))
   , m_nrespa(1)
{
   this->kickoff();
}

void VerletIntegrator::printDetail(FILE* o)
{
   m_thermo->printDetail(o);
   m_baro->printDetail(o);
   printBasic(o);
}

void VerletIntegrator::dynamic(int istep, time_prec dt)
{
   int vers0 = rc_flag & calc::vmask;
   int vers1 = vers0;

   bool save = this->ifSave(istep);
   bool applyBaro = m_baro->ifApply(istep);
   bool mcbaro = false;
   if (applyBaro)
      if (m_baro->getBarostatEnum() == BarostatEnum::MonteCarlo)
         mcbaro = true;
   // toggle off the calc::virial bit if Monte Carlo Barostat is in use
   if (mcbaro)
      vers1 &= ~calc::virial;
   // toggle off the calc::energy bit if neither save nor mcbaro
   if (not save and not mcbaro)
      vers1 &= ~calc::energy;

   time_prec dt2 = 0.5 * dt;
   time_prec dta = dt / m_nrespa;
   time_prec dta2 = 0.5 * dta;

   if (applyBaro)
      m_baro->control1(dt);
   m_thermo->control1(dt);

   if (m_nrespa == 1)
      updateVelocity(dt2);
   else
      updateVelocity2(dta2, dt2);
   this->rattleSave();

   if (applyBaro)
      m_baro->control3(dt);

   if (m_nrespa == 1) {
      updatePosition(dt);
      this->rattle(dt);
      copy_pos_to_xyz(true);
      energy(vers1);
   } else {
      virial_prec vir_fast[9] = {0};
      virial_prec vir_f[9];
      energy_prec esum_f;

      for (int ifast = 1; ifast < m_nrespa; ++ifast) {
         updatePosition(dta);
         copy_pos_to_xyz(false);
         energy(vers1, RESPA_FAST, respa_tsconfig());
         updateVelocity(dta);
         copy_virial(vers1, vir_f);
         if (vers1 & calc::virial)
            for (int i = 0; i < 9; ++i)
               vir_fast[i] += vir_f[i];
      }

      // fast force
      energy(vers1, RESPA_FAST, respa_tsconfig());
      darray::copy(g::q0, n, gx1, gx);
      darray::copy(g::q0, n, gy1, gy);
      darray::copy(g::q0, n, gz1, gz);
      copy_energy(vers1, &esum_f);
      copy_virial(vers1, vir_f);
      if (vers1 & calc::virial)
         for (int i = 0; i < 9; ++i)
            vir_fast[i] += vir_f[i];

      // slow force
      energy(vers1, RESPA_SLOW, respa_tsconfig());
      darray::copy(g::q0, n, gx2, gx);
      darray::copy(g::q0, n, gy2, gy);
      darray::copy(g::q0, n, gz2, gz);
      if (vers1 & calc::energy)
         esum += esum_f;
      if (vers1 & calc::virial)
         for (int i = 0; i < 9; ++i)
            vir[i] += vir_fast[i] / m_nrespa;
   }

   if (applyBaro)
      m_baro->control4(dt);

   if (m_nrespa == 1)
      updateVelocity(dt2);
   else
      updateVelocity2(dta2, dt2);
   this->rattle2(dt, vers1 & calc::virial);

   m_thermo->control2(dt, save);
   if (applyBaro)
      m_baro->control2(dt);
}

void VerletIntegrator::updateVelocity2(time_prec tfast, time_prec tslow)
{
   propagate_velocity2(tfast, gx1, gy1, gz1, tslow, gx2, gy2, gz2);
}
}
