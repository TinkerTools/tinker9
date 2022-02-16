#include "itgbLogVIso.h"
#include "box.h"
#include "mathfunc_ou.h"
#include "mdegv.h"
#include "nose.h"
#include "random.h"
#include "tool/io_print.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void LogVIsoBarostat::control_1_2(time_prec dt)
{
   time_prec dt2 = dt * 0.5;
   double vol0 = volbox();
   double b = 0;
   if (m_langevin)
      b = 2 * m_fric * units::gasconst * bath::kelvin / qbar;

   double dim = 3.0;
   double al = 1.0 + dim / m_dofP;
   double tr_vir = vir[0] + vir[4] + vir[8];

   double gbar =
      2 * eksum * al - tr_vir - dim * vol0 * bath::atmsph / units::prescon;
   gbar /= qbar;
   if (m_langevin)
      vbar = ornstein_uhlenbeck_process(dt2, vbar, m_fric, gbar, b, m_rdn);
   else
      vbar += gbar * dt2;
}

LogVIsoBarostat::LogVIsoBarostat(double fric)
   : BasicBarostat()
   , m_fric(fric)
   , m_rdn()
   , m_langevin(fric != 0.0)
{
   if (m_langevin)
      m_rdn = normal<double>();

   m_dofP = mdstuf::nfree;

   double kt = units::gasconst * bath::kelvin;
   qbar = kt * bath::taupres * bath::taupres * m_dofP;
   vbar = 0;
}

BarostatEnum LogVIsoBarostat::getBarostatEnum() const
{
   return BarostatEnum::m_LogVIso;
}

void LogVIsoBarostat::printDetail(FILE* o)
{
   print(o, " VBar Mass        : %12.2lf\n", qbar);
   printBasic(o);
}

void LogVIsoBarostat::control1(time_prec dt)
{
   if (not applyBaro)
      return;

   control_1_2(dt);
}

void LogVIsoBarostat::control2(time_prec dt)
{
   if (not applyBaro)
      return;

   m_rdn = normal<double>();
   control_1_2(dt);
}

void LogVIsoBarostat::control3(time_prec dt)
{
   if (not applyBaro)
      return;

   double vt = vbar * dt;
   double eterm2 = std::exp(vt);
   lvec1 *= eterm2;
   lvec2 *= eterm2;
   lvec3 *= eterm2;
   set_default_recip_box();
}
}
