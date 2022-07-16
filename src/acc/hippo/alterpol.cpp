//#include "ff/amoebamod.h"
#include "ff/atom.h"
#include "ff/hippomod.h"
#include "ff/image.h"
#include "ff/nblist.h"
//#include "ff/pme.h"
#include "ff/switch.h"
#include "seq/pair_polar_chgpen.h"
//#include "tool/gpucard.h"
#include <array>
#include <fstream>

namespace tinker {
void alterpol()
{
   real cut = switchCut(Switch::REPULS);
   real off = switchOff(Switch::REPULS);

   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;   // question
   const auto* mlst = mlist_unit.deviceptr(); // question

   size_t bufsize = bufferSize(); // question
   PairPolarGrad pgrad;           // question

   // initialize polscale and polinv
   for (int i = 0; i < n; ++i) {
      polscale[i][1][1] = 1.f;
      polscale[i][1][2] = 0.f;
      polscale[i][1][3] = 0.f;
      polscale[i][2][1] = 0.f;
      polscale[i][2][2] = 1.f;
      polscale[i][2][3] = 0.f;
      polscale[i][3][1] = 0.f;
      polscale[i][3][2] = 0.f;
      polscale[i][3][3] = 1.f;
      polinv[i][1][1] = 1.f;
      polinv[i][1][2] = 0.f;
      polinv[i][1][3] = 0.f;
      polinv[i][2][1] = 0.f;
      polinv[i][2][2] = 1.f;
      polinv[i][2][3] = 0.f;
      polinv[i][3][1] = 0.f;
      polinv[i][3][2] = 0.f;
      polinv[i][3][3] = 1.f;
   }

   // find variable polarizability scale matrix at each site
   for (int i = 0; i < n; ++i) {
      real springi = kpep[i];
      real sizi = prepep[i];
      real alphai = dmppep[i];
      int epli = lpep[i];

      int nmlsti = mlst->nlst[i]; // question
      int base = i * maxnlst;     // question
      for (int kk = 0; kk < nmlsti; ++kk) {
         int offset = kk & (bufsize - 1); // question (copied from epolar.cpp)
         int k = mlst->lst[base + kk];    // question
         int eplk = lpep[k];
         if (epli || eplk) {
            real xr = x[k] - x[i];
            real yr = y[k] - y[i];
            real zr = z[k] - z[i];

            zero(pgrad); // question
            real r2 = image2(xr, yr, zr);
            if (r2 <= off2) {
               real r = REAL_SQRT(r2);
               real sizk = prepep[k];
               real alphak = dmppep[k];
               // call dampexpl
            }
         }
      }
   }
}
}