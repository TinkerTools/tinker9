#include "ff/amoebamod.h"
#include "ff/atom.h"
#include "ff/hippomod.h"
#include "ff/image.h"
#include "ff/nblist.h"
//#include "ff/pme.h"
#include "ff/switch.h"
#include "seq/pair_alterpol.h"
#include <tinker/routines.h>
// #include "seq/pair_polar_chgpen.h"
#include "tool/gpucard.h"
// #include <array>
// #include <fstream>

namespace tinker {
void alterpol(real (*polscale)[3][3], real (*polinv)[3][3])
{
   real cut = switchCut(Switch::REPULS);
   real off = switchOff(Switch::REPULS);

   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   // initialize polscale and polinv
   #pragma acc parallel loop independent async\
               deviceptr(polscale)
   for (int i = 0; i < n; ++i) {
      polscale[i][0][0] = 1.f;
      polscale[i][0][1] = 0.f;
      polscale[i][0][2] = 0.f;
      polscale[i][1][0] = 0.f;
      polscale[i][1][1] = 1.f;
      polscale[i][1][2] = 0.f;
      polscale[i][2][0] = 0.f;
      polscale[i][2][1] = 0.f;
      polscale[i][2][2] = 1.f;
   }

   // find variable polarizability scale matrix at each site
   MAYBE_UNUSED int GRID_DIM = gpuGridSize(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(x,y,z,kpep,prepep,dmppep,lpep,mlst,polscale)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      real springi = kpep[i];
      real sizi = prepep[i];
      real alphai = dmppep[i];
      int epli = lpep[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];

      int nmlsti = mlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop vector independent
      for (int kk = 0; kk < nmlsti; ++kk) {
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;
         real r2 = image2(xr, yr, zr);
         int eplk = lpep[k];
         bool incl = (epli || eplk);
         if (r2 <= off2 and incl) {
            real r = REAL_SQRT(r2);
            real springk = kpep[k];
            real sizk = prepep[k];
            real alphak = dmppep[k];
            real ks2i[3][3], ks2k[3][3];
            pair_alterpol(scrtyp, r, r2, 1, cut, off, xr, yr, zr, springi, sizi, alphai, springk,
               sizk, alphak, ks2i, ks2k);
            #pragma acc loop seq
            for (int l = 0; l < 3; ++l) {
               #pragma acc loop seq
               for (int m = 0; m < 3; ++m) {
                  polscale[i][m][l] = polscale[i][m][l] + ks2i[m][l];
                  polscale[k][m][l] = polscale[k][m][l] + ks2k[m][l];
               }
            }
         }
      }
   }

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,kpep,prepep,dmppep,lpep,mlst,mdwexclude,mdwexclude_scale,polscale)
   for (int ii = 0; ii < nmdwexclude; ++ii) {
      int i = mdwexclude[ii][0];
      int k = mdwexclude[ii][1];
      real dscale = mdwexclude_scale[ii][1] - 1;

      real springi = kpep[i];
      real sizi = prepep[i];
      real alphai = dmppep[i];
      int epli = lpep[i];

      real springk = kpep[k];
      real sizk = prepep[k];
      real alphak = dmppep[k];
      int eplk = lpep[k];

      real xr = x[k] - x[i];
      real yr = y[k] - y[i];
      real zr = z[k] - z[i];

      real r2 = image2(xr, yr, zr);
      bool incl1 = dscale != 0;
      bool incl2 = (epli || eplk);
      if (r2 <= off2 and incl1 and incl2) {
         real r = REAL_SQRT(r2);
         real ks2i[3][3], ks2k[3][3];
         pair_alterpol(scrtyp, r, r2, 1, cut, off, xr, yr, zr, springi, sizi, alphai, springk, sizk,
            alphak, ks2i, ks2k);
         #pragma acc loop seq
         for (int l = 0; l < 3; ++l) {
            #pragma acc loop seq
            for (int m = 0; m < 3; ++m) {
               polscale[i][m][l] = polscale[i][m][l] + ks2i[m][l];
               polscale[k][m][l] = polscale[k][m][l] + ks2k[m][l];
            }
         }
      }
   }
   
   // invert
   #pragma acc parallel loop independent async\
               deviceptr(polscale,polinv)
   for (int i = 0; i < n; ++i) {
      real tmp[3][3];
      #pragma acc loop seq
      for (int j = 0; j < 3; ++j) {
         #pragma acc loop seq
         for (int k = 0; k < 3; ++k) {
            tmp[j][k] = polscale[i][j][k];
         }
      }
      real det;
      det = tmp[0][0] * (tmp[1][1] * tmp[2][2] - tmp[1][2] * tmp[2][1]) -
         tmp[1][0] * (tmp[0][1] * tmp[2][2] - tmp[2][1] * tmp[0][2]) +
         tmp[2][0] * (tmp[0][1] * tmp[1][2] - tmp[1][1] * tmp[0][2]);
      polinv[i][0][0] = (tmp[1][1] * tmp[2][2] - tmp[1][2] * tmp[2][1]) / det;
      polinv[i][1][0] = (tmp[2][0] * tmp[1][2] - tmp[1][0] * tmp[2][2]) / det;
      polinv[i][2][0] = (tmp[1][0] * tmp[2][1] - tmp[2][0] * tmp[1][1]) / det;
      polinv[i][0][1] = (tmp[2][1] * tmp[0][2] - tmp[0][1] * tmp[2][2]) / det;
      polinv[i][1][1] = (tmp[0][0] * tmp[2][2] - tmp[2][0] * tmp[0][2]) / det;
      polinv[i][2][1] = (tmp[0][1] * tmp[2][0] - tmp[0][0] * tmp[2][1]) / det;
      polinv[i][0][2] = (tmp[0][1] * tmp[1][2] - tmp[0][2] * tmp[1][1]) / det;
      polinv[i][1][2] = (tmp[0][2] * tmp[1][0] - tmp[0][0] * tmp[1][2]) / det;
      polinv[i][2][2] = (tmp[0][0] * tmp[1][1] - tmp[0][1] * tmp[1][0]) / det;
   }
}

void dexpol() {

}
}
