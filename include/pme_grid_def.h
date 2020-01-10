#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
enum
{
   PCHG_GRID = 1,
   MPOLE_GRID,
   UIND_GRID,
   UIND_GRID_FPHI2,
   DISP_GRID
};
TINKER_NAMESPACE_END


using namespace TINKER_NAMESPACE;
extern "C"
{
   struct PCHG
   {
      static constexpr int N = PCHG_GRID;
   };
   struct MPOLE
   {
      static constexpr int N = MPOLE_GRID;
   };
   struct UIND
   {
      static constexpr int N = UIND_GRID;
   };
   struct UIND2
   {
      static constexpr int N = UIND_GRID_FPHI2;
   };
   struct DISP
   {
      static constexpr int N = DISP_GRID;
   };
}
