#include "add.h"
#include "launch.h"
#include "lpiston.h"
#include "mdegv.h"
#include "mdpq.h"
#include "tool/darray.h"


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
   const mass_prec* restrict molmass,

   const pos_prec* restrict ratcom_x, const pos_prec* restrict ratcom_y,
   const pos_prec* restrict ratcom_z)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;


   for (int im = ithread; im < nmol; im += stride) {
      double mvir = 0;
      double igx, igy, igz;
      double mgx = 0, mgy = 0, mgz = 0;
      vel_prec px = 0, py = 0, pz = 0;
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
         mvir -= (igx * xpos[k] + igy * ypos[k] + igz * zpos[k]);


         mgx += igx;
         mgy += igy;
         mgz += igz;


         mass_prec massk = mass[k];
         px += massk * vx[k];
         py += massk * vy[k];
         pz += massk * vz[k];
      }
      mvir += mgx * ratcom_x[im] + mgy * ratcom_y[im] + mgz * ratcom_z[im];
      auto mmassinv = 1.0 / molmass[im];
      auto mv2 = (px * px + py * py + pz * pz) * mmassinv;


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
              rattle_dmol.molmass, ratcom_x, ratcom_y, ratcom_z);


   val = reduce_sum(lp_molpres_buf, bufsize, g::q0);


   virial_prec atomic_vir = vir[0] + vir[4] + vir[8];
   val += atomic_vir;
}
}
