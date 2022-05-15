#include "ff/amoebamod.h"
#include "seq/launch.h"

namespace tinker {
__global__
static void diagPrecond_cu1(int n, const real* restrict polarity, //
   const real (*restrict rsd)[3], real (*restrict zrsd)[3],       //
   const real (*restrict rsdp)[3], real (*restrict zrsdp)[3])
{
   if (rsdp) {
      for (int i = ITHREAD; i < n; i += STRIDE) {
         real poli = polarity[i];
         zrsd[i][0] = poli * rsd[i][0];
         zrsd[i][1] = poli * rsd[i][1];
         zrsd[i][2] = poli * rsd[i][2];
         zrsdp[i][0] = poli * rsdp[i][0];
         zrsdp[i][1] = poli * rsdp[i][1];
         zrsdp[i][2] = poli * rsdp[i][2];
      }
   } else {
      for (int i = ITHREAD; i < n; i += STRIDE) {
         real poli = polarity[i];
         zrsd[i][0] = poli * rsd[i][0];
         zrsd[i][1] = poli * rsd[i][1];
         zrsd[i][2] = poli * rsd[i][2];
      }
   }
}

void diagPrecond_cu(const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
   launch_k1s(g::s0, n, diagPrecond_cu1, //
      n, polarity,                       //
      rsd, zrsd, rsdp, zrsdp);
}

void diagPrecond2_cu(const real (*rsd)[3], real (*zrsd)[3])
{
   launch_k1s(g::s0, n, diagPrecond_cu1, //
      n, polarity,                       //
      rsd, zrsd, nullptr, nullptr);
}
}
