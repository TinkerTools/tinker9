#pragma once

namespace tinker {
enum class ebond_t
{
   harmonic,
   morse
};

enum class eangle_t : int
{
   in_plane,
   harmonic,
   linear,
   fourier
};

enum class eopbend_t
{
   w_d_c,
   allinger
};
}

extern "C"
{
   // Bond terms.
   struct HARMONIC
   {
      int foo;
   };

   struct MORSE
   {
      int foo;
   };

   // Opbend terms.
   struct WDC
   {
      int foo;
   };

   struct ALLINGER
   {
      int foo;
   };
}
