#include "itgiRespa.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
void RespaIntegrator::kickoff()
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

   energy(calc::grad | calc::virial);
}

RespaIntegrator::~RespaIntegrator()
{
   darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);
}

RespaIntegrator::RespaIntegrator(ThermostatEnum te, BarostatEnum be)
   : VerletIntegrator(te, be)
{
   m_nrespa = mdstuf::nrespa;
   assert(m_nrespa > 1);
   darray::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);
   this->kickoff();
}
}
