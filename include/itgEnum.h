#pragma once

namespace tinker {
enum class PropagatorEnum
{
   Respa,
   Verlet,

   m_LogV,
};

class BasicPropagator;
BasicPropagator* create(PropagatorEnum pe);

enum class ThermostatEnum
{
   Null,
   Andersen,  // empty
   Berendsen, // empty
   Bussi,
   Nhc,

   m_LeapFrogLP,
   m_LP2022,
   m_Nhc1996,
   m_Nhc2006,
};

class BasicThermostat;
BasicThermostat* create(ThermostatEnum);

enum class BarostatEnum
{
   Null,
   Berendsen,
   Bussi, // empty
   LP2022,
   MonteCarlo,
   Nhc2006,

   m_LeapFrogLP,
   m_LogVIso,
   m_LogVAniso,
   m_Nhc1996,
};

class BasicBarostat;
BasicBarostat* create(BarostatEnum);

enum class IntegratorEnum
{
   Beeman,
   Respa,
   Verlet,

   LeapFrogLP,
   LP2022,
   Nhc1996,
   Nhc2006, // J. Phys. A, 39 5629 (2006)
            // https://doi.org/10.1088/0305-4470/39/19/S18
};

class IntegratorStaticData
{
protected:
   static int nrespa;
   static bool applyBaro;
   static double dofP;
   static int anisoArrayLength;
   static const int anisoArray[6][2];
   static constexpr int OrthoOrOct = 3;
   static constexpr int Mono = 4;
   static constexpr int Tri = 6;
};

void __PlaceHolderMessage(const char* s);
}
