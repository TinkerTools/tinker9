#include "itgbMonteCarlo.h"
#include "itgEnum.h"
#include "mdegv.h"
#include "mdpq.h"
#include "mdpt.h"
#include "random.h"
#include "tinker_rt.h"
#include "tool/darray.h"
#include <tinker/detail/bath.hh>

namespace tinker {
MonteCarloBarostat::~MonteCarloBarostat()
{
   darray::deallocate(x_pmonte, y_pmonte, z_pmonte);
}

MonteCarloBarostat::MonteCarloBarostat()
   : BasicBarostat(BarostatEnum::MonteCarlo)
{
   darray::allocate(n, &x_pmonte, &y_pmonte, &z_pmonte);
   get_kv("VOLUME-TRIAL", m_nbaro, 1);
}

void MonteCarloBarostat::control4(time_prec)
{
   T_prec temp = bath::kelvin;
   if (not bath::isothermal)
      kinetic(temp);
   monte_carlo_barostat(esum, temp);
}

bool MonteCarloBarostat::ifApply(int)
{
   double rdm = random<double>();
   if (rdm < 1.0 / m_nbaro)
      return true;
   else
      return false;
}
}
