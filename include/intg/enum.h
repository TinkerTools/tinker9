#pragma once

namespace tinker {
enum class ThermostatEnum
{
   Null,
   Andersen,
   Berendsen,
   Bussi,
   LeapFrogLP,
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
   MonteCarlo,
   Nhc96,
};

class BasicBarostat;
BasicBarostat* create(BarostatEnum);

enum class IntegratorEnum
{
   Beeman,
   LangevinNpt,
   LeapFrogLP,
   Nhc06, // J. Phys. A, 39 5629 (2006)
          // https://doi.org/10.1088/0305-4470/39/19/S18
   Nhc96,
   Verlet,
   VerletRespa,
};
}
