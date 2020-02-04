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
}
