#pragma once
#include "ff/precision.h"
#include "tool/rcman.h"

extern "C"
{
   class PCHG
   {
      int foo;
   };

   class MPOLE
   {
      int foo;
   };

   class UIND
   {
      int foo;
   };

   class UIND2
   {
      int foo;
   };

   class DISP
   {
      int foo;
   };

   class EWALD
   {
      int foo;
   };

   class DEWALD
   {
      int foo;
   };

   class NON_EWALD
   {
      int foo;
   };

   class NON_EWALD_TAPER
   {
      int foo;
   };

   class GORDON1
   {
      int foo;
   };

   class GORDON2
   {
      int foo;
   };
}

namespace tinker {
/// \addtogroup ff
/// \{

bool useEwald();
void elecData(RcOp);

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

TINKER_EXTERN real electric;
TINKER_EXTERN real dielec;
TINKER_EXTERN real elam;

/// \}
}
