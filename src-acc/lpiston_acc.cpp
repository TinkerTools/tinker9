#include "add.h"
#include "lpiston.h"
#include "mdegv.h"
#include "mdpq.h"
#include "rattle.h"


namespace tinker {
void propagate_pos_raxbv_acc(

   pos_prec* r1, pos_prec* r2, pos_prec* r3,

   pos_prec a, pos_prec* x1, pos_prec* x2, pos_prec* x3,

   pos_prec b, pos_prec* y1, pos_prec* y2, pos_prec* y3)
{
   #pragma acc parallel loop independent async\
               deviceptr(r1,r2,r3,x1,x2,x3,y1,y2,y3)
   for (int i = 0; i < n; ++i) {
      r1[i] = r1[i] + a * x1[i] + b * y1[i];
      r2[i] = r2[i] + a * x2[i] + b * y2[i];
      r3[i] = r3[i] + a * x3[i] + b * y3[i];
   }
}


void lp_propagate_mol_vel_acc(vel_prec scal)
{
   auto* molec = rattle_dmol.molecule;
   #pragma acc parallel loop independent async\
               deviceptr(vx,vy,vz,ratcom_vx,ratcom_vy,ratcom_vz,molec)
   for (int i = 0; i < n; ++i) {
      int im = molec[i];
      vx[i] = vx[i] + scal * ratcom_vx[im];
      vy[i] = vy[i] + scal * ratcom_vy[im];
      vz[i] = vz[i] + scal * ratcom_vz[im];
   }
}


void lp_mol_virial_acc()
{
   int nmol = rattle_dmol.nmol;
   auto* molmass = rattle_dmol.molmass;
   auto* imol = rattle_dmol.imol;
   auto* kmol = rattle_dmol.kmol;

   double mvxx = 0, mvyy = 0, mvzz = 0, mvxy = 0, mvxz = 0, mvyz = 0;
   #pragma acc parallel loop independent async\
               copy(mvxx,mvyy,mvzz,mvxy,mvxz,mvyz)\
               reduction(+:mvxx,mvyy,mvzz,mvxy,mvxz,mvyz)\
               deviceptr(molmass,imol,kmol,mass,xpos,ypos,zpos,gx,gy,gz)
   for (int im = 0; im < nmol; ++im) {
      double vxx = 0, vyy = 0, vzz = 0, vxy = 0, vxz = 0, vyz = 0;
      double igx, igy, igz;             // atomic gradients
      pos_prec irx, iry, irz;           // atomic positions
      double mgx = 0, mgy = 0, mgz = 0; // molecular gradients
      pos_prec rx = 0, ry = 0, rz = 0;  // molecular positions
      int start = imol[im][0];
      int end = imol[im][1];
      #pragma acc loop seq
      for (int i = start; i < end; ++i) {
         int k = kmol[i];
#if TINKER_DETERMINISTIC_FORCE
         igx = to_flt_acc<double>(gx[k]);
         igy = to_flt_acc<double>(gy[k]);
         igz = to_flt_acc<double>(gz[k]);
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
      mvxx += vxx;
      mvyy += vyy;
      mvzz += vzz;
      mvxy += vxy;
      mvxz += vxz;
      mvyz += vyz;
   }
   #pragma acc wait

   lp_vir[0] = mvxx + vir[0];
   lp_vir[1] = mvxy + vir[1];
   lp_vir[2] = mvxz + vir[2];
   lp_vir[3] = mvxy + vir[3];
   lp_vir[4] = mvyy + vir[4];
   lp_vir[5] = mvyz + vir[5];
   lp_vir[6] = mvxz + vir[6];
   lp_vir[7] = mvyz + vir[7];
   lp_vir[8] = mvzz + vir[8];
}


void lp_center_of_mass_acc(const pos_prec* ax, const pos_prec* ay,
                           const pos_prec* az, pos_prec* mx, pos_prec* my,
                           pos_prec* mz)
{
   const int nmol = rattle_dmol.nmol;
   const auto* imol = rattle_dmol.imol;
   const auto* kmol = rattle_dmol.kmol;
   const auto* mfrac = ratcom_massfrac;
   #pragma acc parallel loop independent async\
               deviceptr(ax,ay,az,mx,my,mz,mfrac,imol,kmol)
   for (int im = 0; im < nmol; ++im) {
      int start = imol[im][0];
      int end = imol[im][1];
      pos_prec tx = 0, ty = 0, tz = 0;
      #pragma acc loop seq
      for (int i = start; i < end; ++i) {
         int k = kmol[i];
         auto frk = mfrac[k];
         tx += frk * ax[k];
         ty += frk * ay[k];
         tz += frk * az[k];
      }
      mx[im] = tx;
      my[im] = ty;
      mz[im] = tz;
   }
}
}
