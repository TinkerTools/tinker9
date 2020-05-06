#pragma once
#include "compare_types.h"
#include "mdcalc.h"


extern "C"
{
   // PME grids.
   struct PCHG
   {};
   struct MPOLE
   {};
   struct UIND
   {};
   struct UIND2
   {};
   struct DISP
   {};


   // Ewald vs. Non-Ewald
   struct EWALD
   {};
   struct NON_EWALD
   {};
   struct NON_EWALD_TAPER
   {}; // Non-EWALD partial charge also uses switching functions.


   // Energy versions.
   struct Eng : public tinker::EnergyVersion<tinker::calc::v0>
   {};
   struct EngGradVir : public tinker::EnergyVersion<tinker::calc::v1>
   {};
   struct EngAlyz : public tinker::EnergyVersion<tinker::calc::v3>
   {};
   struct EngGrad : public tinker::EnergyVersion<tinker::calc::v4>
   {};
   struct Grad : public tinker::EnergyVersion<tinker::calc::v5>
   {};
   struct GradVir : public tinker::EnergyVersion<tinker::calc::v6>
   {};


   // Bond terms.
   struct HARMONIC
   {};
   struct MORSE
   {};


   // Opbend terms.
   struct WDC
   {};
   struct ALLINGER
   {};


   // VDW terms.
   struct LJ
   {};
   struct BUCK
   {};
   struct MM3HB
   {};
   struct HAL
   {};
   struct GAUSS
   {};
}


namespace tinker {
namespace calc {
using V0 = Eng;
using V1 = EngGradVir;
using V3 = EngAlyz;
using V4 = EngGrad;
using V5 = Grad;
using V6 = GradVir;
}
}
