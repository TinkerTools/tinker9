#pragma once

namespace tinker {
enum class ThermostatEnum
{
   Null,
   Bussi,
   Nhc
};

class BasicThermostat;
BasicThermostat* create(ThermostatEnum);

enum class BarostatEnum
{
   Null,
   MonteCarlo,
   Nhc
};

class BasicBarostat;
BasicBarostat* create(BarostatEnum);

enum class IntegratorEnum
{
   Verlet,
   VerletRespa,
   LeapFrog,
   LeapFrogRespa,
   Nhc06 // J. Phys. A, 39 5629 (2006)
         // https://doi.org/10.1088/0305-4470/39/19/S18
};
}
