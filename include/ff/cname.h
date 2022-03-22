#pragma once
#include "md.h"

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
   struct DEWALD
   {}; // Dispersion PME
   struct NON_EWALD
   {};
   struct NON_EWALD_TAPER
   {}; // Non-EWALD partial charge also uses switching functions.

   // Energy versions.
   struct Eng : public tinker::calc::Vers<tinker::calc::v0>
   {};
   struct EngGradVir : public tinker::calc::Vers<tinker::calc::v1>
   {};
   struct EngAlyz : public tinker::calc::Vers<tinker::calc::v3>
   {};
   struct EngGrad : public tinker::calc::Vers<tinker::calc::v4>
   {};
   struct Grad : public tinker::calc::Vers<tinker::calc::v5>
   {};
   struct GradVir : public tinker::calc::Vers<tinker::calc::v6>
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

   // GORDON1 vs. GORDON2 damping functions
   struct GORDON1
   {};
   struct GORDON2
   {};
}
