#include "ff/atom.h"
#include "ff/image.h"
#include "ff/molecule.h"

namespace tinker {
void copyPosToXyz_acc()
{
   if CONSTEXPR (sizeof(pos_prec) == sizeof(real))
      return;

   #pragma acc parallel loop independent async deviceptr(x,y,z,xpos,ypos,zpos)
   for (int i = 0; i < n; ++i) {
      x[i] = xpos[i];
      y[i] = ypos[i];
      z[i] = zpos[i];
   }
}

void bounds_acc()
{
   auto nmol = molecule.nmol;
   const auto* imol = molecule.imol;
   const auto* kmol = molecule.kmol;
   const auto* molmass = molecule.molmass;

   #pragma acc parallel loop independent async\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(imol,kmol,xpos,ypos,zpos,mass,molmass)
   for (int i = 0; i < nmol; ++i) {
      // locate the center of each molecule
      pos_prec xmid = 0, ymid = 0, zmid = 0;
      int start = imol[i][0];
      int stop = imol[i][1];
      #pragma acc loop seq
      for (int j = start; j < stop; ++j) {
         int k = kmol[j];
         xmid += xpos[k] * mass[k];
         ymid += ypos[k] * mass[k];
         zmid += zpos[k] * mass[k];
      }
      auto weigh = molmass[i];
      xmid /= weigh;
      ymid /= weigh;
      zmid /= weigh;

      // locate the image of the center inside PBC box
      real xc, yc, zc;
      xc = xmid;
      yc = ymid;
      zc = zmid;
      image(xc, yc, zc);

      #pragma acc loop seq
      for (int j = start; j < stop; ++j) {
         int k = kmol[j];
         xpos[k] += xc - xmid;
         ypos[k] += yc - ymid;
         zpos[k] += zc - zmid;
      }
   }
}
}
