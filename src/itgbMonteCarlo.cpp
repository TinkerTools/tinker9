#include "math/random.h"
#include "md/inc.h"
#include "md/integrator.h"
#include "tinker9.h"
#include "tool/darray.h"
#include "tool/io.h"
#include <tinker/detail/bath.hh>

namespace tinker {
MonteCarloBarostat::~MonteCarloBarostat()
{
   darray::deallocate(x_pmonte, y_pmonte, z_pmonte);
}

MonteCarloBarostat::MonteCarloBarostat()
   : BasicBarostat()
{
   darray::allocate(n, &x_pmonte, &y_pmonte, &z_pmonte);
   get_kv("VOLUME-TRIAL", m_nbaro, bath::voltrial);
}

void MonteCarloBarostat::printDetail(FILE* o)
{
   print(o, "\n");
   print(o, " Monte Carlo Barostat\n");
   printBasic(o);
}

BarostatEnum MonteCarloBarostat::getBarostatEnum() const
{
   return BarostatEnum::MonteCarlo;
}

void MonteCarloBarostat::control4(time_prec)
{
   if (not applyBaro)
      return;

   T_prec temp = bath::kelvin;
   if (not bath::isothermal)
      mdKinetic(temp);
   mdMonteCarloBarostat(esum, temp);
}

bool MonteCarloBarostat::ifApply(int)
{
   double rdm = random<double>();
   if (rdm < 1.0 / m_nbaro)
      applyBaro = true;
   else
      applyBaro = false;
   return applyBaro;
}
}
