#include "add.h"
#include "ff/rattle.h"
#include "launch.h"
#include "md/md.h"
#include "tool/darray.h"
#include <tinker/detail/units.hh>

namespace tinker {
__global__
void lp_mol_virial_cu1(virial_buffer restrict hc_vir_buf,

   const double* restrict mass, const pos_prec* restrict xpos, const pos_prec* restrict ypos,
   const pos_prec* restrict zpos, const grad_prec* restrict gx, const grad_prec* restrict gy,
   const grad_prec* restrict gz,

   int nmol, const int (*restrict imol)[2], const int* restrict kmol,
   const double* restrict molmass)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;

   for (int im = ithread; im < nmol; im += stride) {
      double vxx = 0, vyy = 0, vzz = 0, vxy = 0, vxz = 0, vyz = 0;
      double igx, igy, igz;             // atomic gradients
      pos_prec irx, iry, irz;           // atomic positions
      double mgx = 0, mgy = 0, mgz = 0; // molecular gradients
      pos_prec rx = 0, ry = 0, rz = 0;  // molecular positions
      int start = imol[im][0];
      int end = imol[im][1];
      for (int i = start; i < end; ++i) {
         int k = kmol[i];
#if TINKER_DETERMINISTIC_FORCE
         igx = fixedTo<double>(gx[k]);
         igy = fixedTo<double>(gy[k]);
         igz = fixedTo<double>(gz[k]);
#else
         igx = gx[k];
         igy = gy[k];
         igz = gz[k];
#endif
         irx = xpos[k];
         iry = ypos[k];
         irz = zpos[k];
         vxx -= igx * irx;
         vyy -= igy * iry;
         vzz -= igz * irz;
         vxy -= 0.5 * (igx * iry + igy * irx);
         vxz -= 0.5 * (igx * irz + igz * irx);
         vyz -= 0.5 * (igy * irz + igz * iry);

         mgx += igx;
         mgy += igy;
         mgz += igz;
         auto massk = mass[k];
         rx += massk * irx;
         ry += massk * iry;
         rz += massk * irz;
      }
      auto mmassinv = 1 / molmass[im];
      vxx += mgx * rx * mmassinv;
      vyy += mgy * ry * mmassinv;
      vzz += mgz * rz * mmassinv;
      vxy += 0.5 * (mgx * ry + mgy * rx) * mmassinv;
      vxz += 0.5 * (mgx * rz + mgz * rx) * mmassinv;
      vyz += 0.5 * (mgy * rz + mgz * ry) * mmassinv;
      atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, hc_vir_buf, ithread);
   }
}

void lp_mol_virial_cu()
{
   auto bufsize = buffer_size();
   darray::zero(g::q0, bufsize, hc_vir_buf);

   launch_k1b(g::s0, n, lp_mol_virial_cu1,

      hc_vir_buf, mass, xpos, ypos, zpos, gx, gy, gz,

      rattle_dmol.nmol, rattle_dmol.imol, rattle_dmol.kmol, rattle_dmol.molmass);

   virial_reduce(hc_vir, hc_vir_buf);
   for (int iv = 0; iv < 9; ++iv)
      hc_vir[iv] += vir[iv];
}
}
