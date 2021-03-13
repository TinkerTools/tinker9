#include "add.h"
#include "launch.h"
#include "lpiston.h"
#include "mdegv.h"
#include "mdpq.h"
#include "tool/darray.h"
#include <tinker/detail/units.hh>


namespace tinker {
__global__
void lp_molpressure_cu1(
   double* restrict lp_buf, double alpha,

   int n, const mass_prec* restrict mass, const pos_prec* restrict xpos,
   const pos_prec* restrict ypos, const pos_prec* restrict zpos,
   const vel_prec* restrict vx, const vel_prec* restrict vy,
   const vel_prec* restrict vz, const grad_prec* restrict gx,
   const grad_prec* restrict gy, const grad_prec* restrict gz,

   int nmol, const int (*restrict imol)[2], const int* restrict kmol,
   const mass_prec* restrict molmass)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;


   for (int im = ithread; im < nmol; im += stride) {
      double mvir = 0;
      double igx, igy, igz;             // atomic gradients
      pos_prec irx, iry, irz;           // atomic positions
      double mgx = 0, mgy = 0, mgz = 0; // molecular gradients
      pos_prec rx = 0, ry = 0, rz = 0;  // molecular positions
      vel_prec px = 0, py = 0, pz = 0;  // molecular momenta
      int start = imol[im][0];
      int end = imol[im][1];
      for (int i = start; i < end; ++i) {
         int k = kmol[i];
#if TINKER_DETERMINISTIC_FORCE
         igx = to_flt_cu<double>(gx[k]);
         igy = to_flt_cu<double>(gy[k]);
         igz = to_flt_cu<double>(gz[k]);
#else
         igx = gx[k];
         igy = gy[k];
         igz = gz[k];
#endif
         irx = xpos[k];
         iry = ypos[k];
         irz = zpos[k];
         mvir -= (igx * irx + igy * iry + igz * irz); // unit: kcal/mol


         mgx += igx;
         mgy += igy;
         mgz += igz;
         mass_prec massk = mass[k];
         rx += massk * irx;
         ry += massk * iry;
         rz += massk * irz;
         px += massk * vx[k];
         py += massk * vy[k];
         pz += massk * vz[k];
      }
      auto mmassinv = 1.0 / molmass[im];
      mvir += (mgx * rx + mgy * ry + mgz * rz) * mmassinv;
      auto mv2 = (px * px + py * py + pz * pz) * mmassinv / units::ekcal;


      lp_buf[ithread] += (alpha * mv2 - mvir);
   }
}


void lp_molpressure_cu(double alpha, double& val)
{
   auto bufsize = buffer_size();
   darray::zero(g::q0, bufsize, lp_molpres_buf);


   launch_k1b(g::s0, n, lp_molpressure_cu1,

              lp_molpres_buf, alpha,

              n, mass, xpos, ypos, zpos, vx, vy, vz, gx, gy, gz,

              rattle_dmol.nmol, rattle_dmol.imol, rattle_dmol.kmol,
              rattle_dmol.molmass);


   double sum = reduce_sum(lp_molpres_buf, bufsize, g::q0);
   virial_prec atomic_vir = vir[0] + vir[4] + vir[8];
   val = sum - atomic_vir;
}
}
