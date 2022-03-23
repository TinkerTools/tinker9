#pragma once

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
