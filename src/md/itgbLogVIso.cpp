#include "ff/box.h"
#include "ff/energy.h"
#include "math/ou.h"
#include "math/random.h"
#include "md/integrator.h"
#include "md/pt.h"
#include "md/rattle.h"
#include "tool/ioprint.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void IsoBaroDevice::control_1_2(time_prec dt)
{
   time_prec dt2 = dt * 0.5;
   double vol0 = boxVolume();
   double b = 0;
   if (m_langevin)
      b = std::sqrt(2 * m_fric * units::gasconst * bath::kelvin / qbar);

   double dim = 3.0;
   double al = 1.0 + dim / dofP;
   double tr_vir = m_vir[0] + m_vir[4] + m_vir[8];
   double eksu1 = *m_eksum;

   double gbar = 2 * eksu1 * al - tr_vir - dim * vol0 * bath::atmsph / units::prescon;
   gbar /= qbar;
   if (m_langevin)
      vbar = OUProcess(dt2, vbar, m_fric, gbar, b, m_rdn);
   else
      vbar += gbar * dt2;
}

IsoBaroDevice::IsoBaroDevice(double fric)
   : BasicBarostat()
   , m_vir(nullptr)
   , m_eksum(nullptr)
   , m_fric(fric)
   , m_rdn()
   , m_langevin(fric != 0.0)
{
   if (atomic) {
      m_vir = vir;
      m_eksum = &eksum;
   } else {
      m_vir = hc_vir;
      m_eksum = &hc_eksum;
   }

   if (m_langevin)
      m_rdn = normal<double>();

   dofP = mdstuf::nfree;

   double kt = units::gasconst * bath::kelvin;
   qbar = kt * bath::taupres * bath::taupres * dofP;
   vbar = 0;
}

BarostatEnum IsoBaroDevice::getBarostatEnum() const
{
   return BarostatEnum::m_LogVIso;
}

void IsoBaroDevice::printDetail(FILE* o)
{
   print(o, " VBar Mass          %12.2lf\n", qbar);
   printBasic(o);
}

void IsoBaroDevice::control1(time_prec dt)
{
   if (not applyBaro)
      return;

   control_1_2(dt);
}

void IsoBaroDevice::control2(time_prec dt)
{
   if (not applyBaro)
      return;

   m_rdn = normal<double>();
   control_1_2(dt);
}

void IsoBaroDevice::control3(time_prec dt)
{
   if (not applyBaro)
      return;

   double vt = vbar * dt;
   double eterm2 = std::exp(vt);
   lvec1 *= eterm2;
   lvec2 *= eterm2;
   lvec3 *= eterm2;
   boxSetCurrentRecip();
}
}
