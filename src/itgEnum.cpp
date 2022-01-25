#include "itgEnum.h"
#include "itgbBasic.h"
#include "itgbBerendsen.h"
#include "itgbMonteCarlo.h"
#include "itgtBasic.h"
#include "itgtBussi.h"
#include "itgtNhc96.h"
#include "mdegv.h"
#include "mdpq.h"
#include "mdpt.h"

namespace tinker {
static double* localAtomKinetic()
{
   T_prec temp;
   kinetic(temp);
   return &eksum;
}

static void localScaleAtomVelocity(double velsc)
{
   darray::scale(g::q0, n, velsc, vx);
   darray::scale(g::q0, n, velsc, vy);
   darray::scale(g::q0, n, velsc, vz);
}

BasicThermostat* create(ThermostatEnum te)
{
   if (te == ThermostatEnum::Bussi)
      return new BussiThermostat;
   else if (te == ThermostatEnum::Nhc96)
      return new Nhc96Thermostat(5, 5, 3.0 * (n - 1), localAtomKinetic,
                                 localScaleAtomVelocity, std::string("NHC"));
   else
      return new BasicThermostat;
}

BasicBarostat* create(BarostatEnum be)
{
   if (be == BarostatEnum::Berendsen)
      return new BerendsenBarostat;
   else if (be == BarostatEnum::MonteCarlo)
      return new MonteCarloBarostat;
   else
      return new BasicBarostat;
}
}
