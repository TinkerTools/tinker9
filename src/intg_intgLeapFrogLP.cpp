#include "energy.h"
#include "intg/intgLeapFrogLP.h"
#include "lf_lpiston.h"
#include "mdcalc.h"
#include "mdpq.h"
#include "tool/darray.h"

namespace tinker {
LeapFrogLPIntegrator::~LeapFrogLPIntegrator()
{
   darray::deallocate(leapfrog_x, leapfrog_y, leapfrog_z);
   darray::deallocate(leapfrog_vx, leapfrog_vy, leapfrog_vz, leapfrog_vxold,
                      leapfrog_vyold, leapfrog_vzold);
}

LeapFrogLPIntegrator::LeapFrogLPIntegrator()
   : BasicIntegrator()
{
   darray::allocate(n, &leapfrog_x, &leapfrog_y, &leapfrog_z);
   darray::allocate(n, &leapfrog_vx, &leapfrog_vy, &leapfrog_vz,
                    &leapfrog_vxold, &leapfrog_vyold, &leapfrog_vzold);
}

void LeapFrogLPIntegrator::printDetail(FILE*) {}

void LeapFrogLPIntegrator::kickoff()
{
   energy(calc::energy | calc::grad | calc::virial);
}

void LeapFrogLPIntegrator::dynamic(int istep, time_prec dt)
{
   lf_lpiston_npt(istep, dt);
}
}
