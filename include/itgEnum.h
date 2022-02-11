#pragma once

namespace tinker {
enum class ThermostatEnum
{
   Null,
   Andersen,
   Berendsen,
   Bussi,
   LeapFrogLP,
   LP22,
   Nhc06,
   Nhc96,
};

class BasicThermostat;
BasicThermostat* create(ThermostatEnum);

enum class BarostatEnum
{
   Null,
   Berendsen,
   Bussi,
   Langevin,
   LeapFrogLP,
   LP22,
   MonteCarlo,
   Nhc06,
   Nhc96,
};

class BasicBarostat;
BasicBarostat* create(BarostatEnum);

enum class IntegratorEnum
{
   Beeman,
   LangevinNpt,
   LeapFrogLP,
   LP22,
   Nhc06, // J. Phys. A, 39 5629 (2006)
          // https://doi.org/10.1088/0305-4470/39/19/S18
   Nhc96,
   Respa,
   Verlet,
};
}
