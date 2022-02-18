#include "itgiNhc96.h"
#include "energy.h"
#include "mdcalc.h"
#include "nose.h"
#include "tool/error.h"
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
const char* Nhc96Integrator::name() const
{
   return "Molecular Dynamics Trajectory via Nose-Hoover NPT Algorithm";
}

void Nhc96Integrator::kickoff()
{
   double ekt = units::gasconst * bath::kelvin;
   vbar = 0;
   qbar = (mdstuf::nfree + 1) * ekt * bath::taupres * bath::taupres;
   gbar = 0;
   for (int i = 0; i < maxnose; ++i) {
      vnh[i] = 0;
      qnh[i] = ekt * bath::tautemp * bath::tautemp;
      gnh[i] = 0;
   }
   qnh[0] *= mdstuf::nfree;
   energy(calc::grad | calc::virial);
}

Nhc96Integrator::Nhc96Integrator()
   : BasicIntegrator()
{
   if (m_prop->useRattle())
      TINKER_THROW("Constraints under NH-NPT require the ROLL algorithm.");

   this->kickoff();
}

void Nhc96Integrator::dynamic(int istep, time_prec dt)
{
   nhc_npt(istep, dt);
}
}
