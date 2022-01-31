#include "itgiNhc06.h"
#include "box.h"
#include "energy.h"
#include "itgEnum.h"
#include "mathfunc_sinhc.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "nose.h"
#include "tinker_rt.h"
#include "tool/error.h"
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
double* Nhc06Thermostat::vbarKinetic()
{
   static double ekvbar;
   ekvbar = 0.5 * qbar * vbar * vbar;
   return &ekvbar;
}

void Nhc06Thermostat::scaleVbarVelocity(double scale)
{
   vbar *= scale;
}

Nhc06Thermostat::Nhc06Thermostat()
   : BasicThermostat(ThermostatEnum::Nhc06)
{
   double dofT = mdstuf::nfree;
   int nhclen = 4;
   int nc = 4;

   // tpart
   m_tpart =
      new Nhc96Thermostat(nhclen, nc, dofT, Nhc96Thermostat::atomicKinetic,
                          Nhc96Thermostat::scaleAtomicVelocity, "NHC");

   // tbaro
   m_tbaro = new Nhc96Thermostat(nhclen, nc, 1.0, Nhc06Thermostat::vbarKinetic,
                                 Nhc06Thermostat::scaleVbarVelocity, "NHCB");
}

Nhc06Thermostat::~Nhc06Thermostat()
{
   delete m_tpart;
   delete m_tbaro;
}

void Nhc06Thermostat::printDetail(FILE* o)
{
   m_tpart->printDetail(o);
   m_tbaro->printDetail(o);
}

void Nhc06Thermostat::control1b(time_prec dt, bool applyBaro)
{
   if (applyBaro)
      m_tbaro->control1(dt);
   m_tpart->control1(dt);
}

void Nhc06Thermostat::control2b(time_prec dt, bool save, bool applyBaro)
{
   m_tpart->control2(dt, save);
   if (applyBaro)
      m_tbaro->control2(dt, false);
}

//====================================================================//

void Nhc06Barostat::control12Impl(time_prec dt)
{
   double dim = 3.0;
   double al = 1.0 + dim / m_dofP;
   double vol0 = volbox();
   double tr_vir = vir[0] + vir[4] + vir[8];
   double gbar =
      al * 2 * eksum - tr_vir - dim * vol0 * bath::atmsph / units::prescon;
   gbar /= qbar;
   vbar += gbar * dt * 0.5;
}

Nhc06Barostat::Nhc06Barostat()
   : BasicBarostat(BarostatEnum::Nhc06)
{
   m_dofP = mdstuf::nfree;

   double kt = units::gasconst * bath::kelvin;
   qbar = kt * bath::taupres * bath::taupres * m_dofP;
   vbar = 0;
}

void Nhc06Barostat::printDetail(FILE* o)
{
   print(o, " VBar Mass        : %12.2lf\n", qbar);
   printBasic(o);
}

double Nhc06Barostat::dof() const
{
   return m_dofP;
}

void Nhc06Barostat::control1(time_prec dt)
{
   control12Impl(dt);
}

void Nhc06Barostat::control2(time_prec dt)
{
   control12Impl(dt);
}

void Nhc06Barostat::control3(time_prec dt)
{
   if (not m_apply)
      return;

   double vt = vbar * dt;
   double eterm2 = std::exp(vt);
   lvec1 *= eterm2;
   lvec2 *= eterm2;
   lvec3 *= eterm2;
   set_default_recip_box();
}

//====================================================================//

void Nhc06Integrator::kickoff()
{
   energy(calc::grad | calc::virial);
}

Nhc06Integrator::Nhc06Integrator()
   : BasicIntegrator()
{
   if (m_userat)
      TINKER_THROW("Constraints under NH-NPT require the ROLL algorithm.");

   m_thermo = new Nhc06Thermostat;
   m_baro = new Nhc06Barostat;
   get_kbool("PEDANTIC", m_pedantic, false);

   this->kickoff();
}

Nhc06Integrator::~Nhc06Integrator()
{
   delete m_thermo;
   delete m_baro;
}

void Nhc06Integrator::printDetail(FILE* o)
{
   m_thermo->printDetail(o);
   m_baro->printDetail(o);
   printBasic(o);
}

void Nhc06Integrator::dynamic(int istep, time_prec dt)
{
   int vers0 = rc_flag & calc::vmask;
   int vers1 = vers0;

   bool save = this->ifSave(istep);
   bool applyBaro = m_baro->ifApply(istep);
   if (not applyBaro)
      vers1 &= ~calc::virial;
   if (not save)
      vers1 &= ~calc::energy;

   time_prec dt2 = dt / 2;

   double al = 0.0;
   if (applyBaro)
      al = 1.0 + 3.0 / m_baro->dof();

   m_thermo->control1b(dt, applyBaro);
   m_baro->control1(dt); // update vbar

   if (m_pedantic) {
      updateVelocityPedantic(dt2, al * vbar);
   } else {
      double scale = std::exp(-al * vbar * dt2);
      Nhc96Thermostat::scaleAtomicVelocity(scale);
      updateVelocity(dt2);
   }

   m_baro->control3(dt); // update the volume

   if (m_pedantic)
      updatePositionPedantic(dt);
   else
      updatePosition(dt);
   copy_pos_to_xyz(true);
   energy(vers1);

   if (m_pedantic) {
      updateVelocityPedantic(dt2, al * vbar);
   } else {
      updateVelocity(dt2);
      double scale = std::exp(-al * vbar * dt2);
      Nhc96Thermostat::scaleAtomicVelocity(scale);
   }

   m_baro->control2(dt); // update vbar
   m_thermo->control2b(dt, save, applyBaro);
}

void Nhc06Integrator::updatePositionPedantic(time_prec t)
{
   double vt2 = vbar * t / 2;
   double expvt2 = std::exp(vt2);
   double coefr = expvt2 * expvt2;
   double coefv = t * expvt2 * sinhc(vt2);
   propagate_pos_axbv(coefr, coefv);
}

void Nhc06Integrator::updateVelocityPedantic(time_prec t, double velbar)
{
   double vt2 = velbar * t / 2;
   double expvt2 = std::exp(-vt2);
   double coefv = expvt2 * expvt2;
   double coeff = t * expvt2 * sinhc(vt2);
   propagate_vel_avbf(coefv, coeff, gx, gy, gz);
}
}
