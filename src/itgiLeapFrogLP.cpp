#include "ff/energy.h"
#include "md/integrator.h"
#include "md/lflpiston.h"
#include "md/md.h"
#include "tool/darray.h"

namespace tinker {
const char* LeapFrogLPIntegrator::name() const
{
   return "Molecular Dynamics Trajectory via Langevin Piston Algorithm";
}

void LeapFrogLPIntegrator::kickoff()
{
   energy(calc::energy | calc::grad | calc::virial);
}

LeapFrogLPIntegrator::LeapFrogLPIntegrator()
   : BasicIntegrator()
{
   darray::allocate(n, &leapfrog_x, &leapfrog_y, &leapfrog_z);
   darray::allocate(n, &leapfrog_vx, &leapfrog_vy, &leapfrog_vz, &leapfrog_vxold, &leapfrog_vyold,
      &leapfrog_vzold);
   this->kickoff();
}

LeapFrogLPIntegrator::~LeapFrogLPIntegrator()
{
   darray::deallocate(leapfrog_x, leapfrog_y, leapfrog_z);
   darray::deallocate(
      leapfrog_vx, leapfrog_vy, leapfrog_vz, leapfrog_vxold, leapfrog_vyold, leapfrog_vzold);
}

void LeapFrogLPIntegrator::dynamic(int istep, time_prec dt)
{
   lf_lpiston_npt(istep, dt);
}
}
