#pragma once
#include "ff/precision.h"
#include "tool/rcman.h"

TINKER_DECL_EXTN("C")
{
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

   struct EWALD
   {
      int foo;
   };

   struct DEWALD
   {
      int foo;
   };

   struct NON_EWALD
   {
      int foo;
   };

   struct NON_EWALD_TAPER
   {
      int foo;
   };

   struct GORDON1
   {
      int foo;
   };

   struct GORDON2
   {
      int foo;
   };
}

namespace tinker {
/// \ingroup ff
bool useEwald();
/// \ingroup ff
void elecData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup ff
TINKER_EXTERN real electric;
/// \ingroup ff
TINKER_EXTERN real dielec;
}
