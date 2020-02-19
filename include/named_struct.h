#pragma once
#include "compare_types.h"
#include "md_calc.h"


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


   // Energy versions.
   struct Eng
      : public TINKER_NAMESPACE::EnergyVersion<TINKER_NAMESPACE::calc::v0>
   {};
   struct EngGradVir
      : public TINKER_NAMESPACE::EnergyVersion<TINKER_NAMESPACE::calc::v1>
   {};
   struct EngAlyz
      : public TINKER_NAMESPACE::EnergyVersion<TINKER_NAMESPACE::calc::v3>
   {};
   struct EngGrad
      : public TINKER_NAMESPACE::EnergyVersion<TINKER_NAMESPACE::calc::v4>
   {};
   struct Grad
      : public TINKER_NAMESPACE::EnergyVersion<TINKER_NAMESPACE::calc::v5>
   {};
   struct GradVir
      : public TINKER_NAMESPACE::EnergyVersion<TINKER_NAMESPACE::calc::v6>
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


TINKER_NAMESPACE_BEGIN
namespace calc {
using V0 = Eng;
using V1 = EngGradVir;
using V3 = EngAlyz;
using V4 = EngGrad;
using V5 = Grad;
using V6 = GradVir;
}
TINKER_NAMESPACE_END
