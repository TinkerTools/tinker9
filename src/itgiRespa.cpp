#include "itgiRespa.h"
#include "energy.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdprec.h"

namespace tinker {
const char* RespaIntegrator::name() const
{
   return "Molecular Dynamics Trajectory via r-RESPA MTS Algorithm";
}

void RespaIntegrator::kickoff()
{
   RespaIntegrator::KickOff();
}

RespaIntegrator::RespaIntegrator(ThermostatEnum te, BarostatEnum be)
   : BasicIntegrator(PropagatorEnum::Respa, te, be)
{
   this->kickoff();
}

void RespaIntegrator::KickOff()
{
   // save fast gradients to gx1 etc.
   energy(calc::grad, RESPA_FAST, respa_tsconfig());
   darray::copy(g::q0, n, gx1, gx);
   darray::copy(g::q0, n, gy1, gy);
   darray::copy(g::q0, n, gz1, gz);

   // save slow gradients to gx2 etc.
   energy(calc::grad, RESPA_SLOW, respa_tsconfig());
   darray::copy(g::q0, n, gx2, gx);
   darray::copy(g::q0, n, gy2, gy);
   darray::copy(g::q0, n, gz2, gz);

   energy((calc::grad | calc::virial) & rc_flag);
}
}
