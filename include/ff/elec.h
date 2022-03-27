#pragma once
#include "precision.h"
#include "tool/rcman.h"

namespace tinker {
bool useEwald();
void elecData(RcOp);
}

extern "C"
{
   // PME grids.
   struct PCHG
   {
      int foo;
   };

   struct MPOLE
   {
      int foo;
   };

   struct UIND
   {
      int foo;
   };

   struct UIND2
   {
      int foo;
   };

   struct DISP
   {
      int foo;
   };

   // Ewald vs. Non-Ewald
   struct EWALD
   {
      int foo;
   };

   struct DEWALD
   {
      int foo;
   }; // Dispersion PME

   struct NON_EWALD
   {
      int foo;
   };

   struct NON_EWALD_TAPER
   {
      int foo;
   }; // Non-EWALD partial charge also uses switching functions.

   // GORDON1 vs. GORDON2 damping functions
   struct GORDON1
   {
      int foo;
   };

   struct GORDON2
   {
      int foo;
   };
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
TINKER_EXTERN real electric;
TINKER_EXTERN real dielec;
}
