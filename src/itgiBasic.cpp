#include "itgiBasic.h"
#include "energy.h"
#include "lpiston.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "tool/io_print.h"

namespace tinker {
void BasicIntegrator::plan(int istep)
{
   int vers0 = rc_flag & calc::vmask;
   vers1 = vers0;

   save = m_prop->ifSave(istep);
   bool mcbaro = false;
   if (m_baro->ifApply(istep))
      if (m_baro->getBarostatEnum() == BarostatEnum::MonteCarlo)
         mcbaro = true;
   // toggle off virial if !applyBaro
   if (not applyBaro)
      vers1 &= ~calc::virial;
   // toggle off virial for MC barostat
   if (mcbaro)
      vers1 &= ~calc::virial;
   // toggle off energy if neither save nor mcbaro
   if (not save and not mcbaro)
      vers1 &= ~calc::energy;
}

BasicIntegrator::BasicIntegrator(PropagatorEnum pe, ThermostatEnum te,
                                 BarostatEnum be)
   : m_prop(create(pe))
   , m_thermo(create(te))
   , m_baro(create(be))
{
   this->plan(0);
}

BasicIntegrator::BasicIntegrator()
   : m_prop(new BasicPropagator)
   , m_thermo(new BasicThermostat)
   , m_baro(new BasicBarostat)
{
   this->plan(0);
}

BasicIntegrator::~BasicIntegrator()
{
   delete m_prop;
   delete m_thermo;
   delete m_baro;
}

void BasicIntegrator::printDetail(FILE* o)
{
   m_thermo->printDetail(o);
   m_baro->printDetail(o);
   print(o, "\n");
   print(o, " NRespa             %12d\n", nrespa);
   print(o, "\n");
   print(o, " %s\n", this->name());
}

void BasicIntegrator::dynamic(int istep, time_prec dt)
{
   time_prec dt2 = dt * 0.5;

   m_baro->control1(dt);
   m_thermo->control1(dt);

   if (nrespa == 1)
      m_prop->updateVelocity1(dt2);
   else
      m_prop->updateVelocityR1(dt2, nrespa);
   m_prop->rattleSave();

   m_baro->control3(dt);

   this->plan(istep);
   if (nrespa == 1) {
      m_prop->updatePosition(dt);
      m_prop->rattle(dt);
      copy_pos_to_xyz(true);
      energy(vers1);
      if (vers1 & calc::virial)
         if (not atomic)
            lp_virial(true);
   } else {
      virial_prec vir_fast[9] = {0};
      energy_prec esum_f;
      time_prec dta = dt / nrespa;

      for (int ifast = 1; ifast < nrespa; ++ifast) {
         m_prop->updatePosition(dta);
         copy_pos_to_xyz(false);
         energy(vers1, RESPA_FAST, respa_tsconfig());
         m_prop->updateVelocityR0(dta);
         if (vers1 & calc::virial) {
            if (atomic) {
               for (int iv = 0; iv < 9; ++iv)
                  vir_fast[iv] += vir[iv];
            } else {
               lp_virial(true);
               for (int iv = 0; iv < 9; ++iv)
                  vir_fast[iv] += lp_vir[iv];
            }
         }
      }
      m_prop->updatePosition(dta);
      m_prop->rattle(dt);
      copy_pos_to_xyz(true);

      // fast force
      energy(vers1, RESPA_FAST, respa_tsconfig());
      darray::copy(g::q0, n, gx1, gx);
      darray::copy(g::q0, n, gy1, gy);
      darray::copy(g::q0, n, gz1, gz);
      copy_energy(vers1, &esum_f);
      if (vers1 & calc::virial) {
         if (atomic) {
            for (int iv = 0; iv < 9; ++iv)
               vir_fast[iv] += vir[iv];
         } else {
            lp_virial(true);
            for (int iv = 0; iv < 9; ++iv)
               vir_fast[iv] += lp_vir[iv];
         }
      }

      // slow force
      energy(vers1, RESPA_SLOW, respa_tsconfig());
      darray::copy(g::q0, n, gx2, gx);
      darray::copy(g::q0, n, gy2, gy);
      darray::copy(g::q0, n, gz2, gz);
      if (vers1 & calc::energy)
         esum += esum_f;
      if (vers1 & calc::virial) {
         if (atomic) {
            for (int iv = 0; iv < 9; ++iv)
               vir[iv] += vir_fast[iv] / nrespa;
         } else {
            lp_virial(true);
            for (int iv = 0; iv < 9; ++iv)
               lp_vir[iv] += vir_fast[iv] / nrespa;
         }
      }
   }

   m_baro->control4(dt);

   if (nrespa == 1)
      m_prop->updateVelocity2(dt2);
   else
      m_prop->updateVelocityR2(dt2, nrespa);
   m_prop->rattle2(dt, vers1 & calc::virial);

   m_thermo->control2(dt, save);
   m_baro->control2(dt);
}
}
