#pragma once
#include "intg/baroMonteCarlo.h"
#include "mdegv.h"
#include "mdpt.h"
#include "random.h"
#include <tinker/detail/bath.hh>

namespace tinker {
MonteCarloBarostat::~MonteCarloBarostat() {}

MonteCarloBarostat::MonteCarloBarostat()
   : BasicBarostat()
{}

void MonteCarloBarostat::control4(time_prec)
{
   T_prec temp = bath::kelvin;
   if (not bath::isothermal)
      kinetic(temp);
   monte_carlo_barostat(esum, temp);
}

bool MonteCarloBarostat::ifApply(int istep)
{
   double rdm = random<double>();
   if (rdm < 1.0 / m_nbaro)
      return true;
   else
      return false;
}
}
