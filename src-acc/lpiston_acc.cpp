#include "add.h"
#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "mdpt.h"
#include "nose.h"
#include "random.h"
#include "rattle.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {


//====================================================================//


void ratcom_kevir_acc(double coef, double atomic_vir, double& val)
{
   int nmol = rattle_dmol.nmol;
   auto* molmass = rattle_dmol.molmass;
   auto* imol = rattle_dmol.imol;
   auto* kmol = rattle_dmol.kmol;

   double half_coeff = 0.5 * coef;

   double sum = 0;
   #pragma acc parallel loop independent async\
               copy(sum) reduction(+:sum)\
               deviceptr(molmass,imol,kmol,ratcom_px,ratcom_py,ratcom_pz,\
                  ratcom_x,ratcom_y,ratcom_z,ratcom_gx,ratcom_gy,ratcom_gz,\
                  xpos,ypos,zpos,gx,gy,gz)
   for (int im = 0; im < nmol; ++im) {
      auto px = ratcom_px[im];
      auto py = ratcom_py[im];
      auto pz = ratcom_pz[im];
      auto mmassinv = 1.0 / molmass[im];
      auto ke_x2 = (px * px + py * py + pz * pz) * mmassinv;

      double igx, igy, igz, mgx, mgy, mgz;
      double mvir = 0.0;
#if TINKER_DETERMINISTIC_FORCE
      mgx = to_flt_acc<double>(ratcom_gx[im]);
      mgy = to_flt_acc<double>(ratcom_gy[im]);
      mgz = to_flt_acc<double>(ratcom_gz[im]);
#else
      mgx = ratcom_gx[im];
      mgy = ratcom_gy[im];
      mgz = ratcom_gz[im];
#endif
      mvir += mgx * ratcom_x[im] + mgy * ratcom_y[im] + mgz * ratcom_z[im];

      int start = imol[im][0];
      int stop = imol[im][1];
      #pragma acc loop seq
      for (int i = start; i < stop; ++i) {
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
         mvir -= (igx * xpos[k] + igy * ypos[k] + igz * zpos[k]);
      }

      sum += half_coeff * ke_x2 - mvir;
   }
   #pragma acc wait

   // sum = coef * kinetic + molecular_vir - atomic_vir
   // return: sum + atomic_vir
   val = sum + atomic_vir;
}


namespace {
void vvlp_hc_pt(time_prec dt, double R, int)
{
   const time_prec dt2 = 0.5 * dt;
   const time_prec dt4 = 0.25 * dt;
   const time_prec dt8 = 0.125 * dt;
   const double D = 3.0;
   const double Nf = mdstuf::nfree;
   const double g = stodyn::friction;
   const double kbt = units::gasconst * bath::kelvin;
   const double gkbt = Nf * kbt;
   double odnf = 1.0;
   if (n > 1)
      odnf = 1.0 + D / (3.0 * (n - 1));
   const double sdbar = std::sqrt(2.0 * kbt * g * dt / qbar) / 2;


   T_prec temp;
   kinetic(temp);


   for (int i = maxnose - 1; i > -1; --i) {
      if (i == 0)
         gnh[i] = (2 * eksum - gkbt) / qnh[i];
      else
         gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];


      if (i == maxnose - 1)
         vnh[i] += gnh[i] * dt4;
      else {
         double exptm = std::exp(-vnh[i + 1] * dt8);
         vnh[i] = (vnh[i] * exptm + gnh[i] * dt4) * exptm;
      }
   }


   gbar = ratcom_kevir_value / qbar;
   vbar += gbar * dt2 - g * vbar * dt2 + sdbar * R;


   for (int i = 0; i < maxnose; ++i) {
      if (i == 0)
         gnh[i] = (2 * eksum - gkbt) / qnh[i];
      else
         gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];


      if (i == maxnose - 1)
         vnh[i] += gnh[i] * dt4;
      else {
         double exptm = std::exp(-vnh[i + 1] * dt8);
         vnh[i] = (vnh[i] * exptm + gnh[i] * dt4) * exptm;
      }
   }
}
}


void vvlp_hc_acc(int istep, time_prec dt)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(istep % inform::iwrite);
   if (!save)
      vers1 &= ~calc::energy;


   const time_prec dt2 = 0.5 * dt;
   const double D = 3.0;
   const double odnf = 1.0 + D / (3 * (n - 1));
   const double R = normal<double>();


   const int nmol = rattle_dmol.nmol;
   const auto* imol = rattle_dmol.imol;
   const auto* kmol = rattle_dmol.kmol;
   const auto* molec = rattle_dmol.molecule;
   const auto* molmass = rattle_dmol.molmass;
   const auto* mfrac = ratcom_massfrac;


   // thermostat 1/2
   vvlp_hc_pt(dt, R, 1);
   const double vnh0 = vnh[0];


   // velocity 1/2
   // iL6 dt/2
   #pragma acc parallel loop independent async\
               deviceptr(vx,vy,vz,mfrac)
   for (int i = 0; i < n; ++i) {
      auto wui = mfrac[i];
      auto coef = std::exp(-dt2 * (vnh0 + odnf * vbar * wui));
      vx[i] *= coef;
      vy[i] *= coef;
      vz[i] *= coef;
   }
   // update Pu
   #pragma acc parallel loop independent async\
               deviceptr(vx,vy,vz,ratcom_px,ratcom_py,ratcom_pz,\
                  mass,imol,kmol)
   for (int im = 0; im < nmol; ++im) {
      int start = imol[im][0];
      int end = imol[im][1];
      vel_prec rpx = 0, rpy = 0, rpz = 0;
      #pragma acc loop seq
      for (int i = start; i < end; ++i) {
         int k = kmol[i];
         auto massk = mass[k];
         rpx += massk * vx[k];
         rpy += massk * vy[k];
         rpz += massk * vz[k];
      }
      ratcom_px[im] = rpx;
      ratcom_py[im] = rpy;
      ratcom_pz[im] = rpz;
   }
   // iL5 dt/2
   #pragma acc parallel loop independent async\
               deviceptr(vx,vy,vz,ratcom_px,ratcom_py,ratcom_pz,\
                  mfrac,molmass,molec)
   for (int i = 0; i < n; ++i) {
      auto im = molec[i];
      auto invMu = 1.0 / molmass[im];
      auto wui = mfrac[i];
      vx[i] -= dt2 * odnf * vbar * (ratcom_px[i] * invMu - wui * vx[i]);
      vy[i] -= dt2 * odnf * vbar * (ratcom_py[i] * invMu - wui * vy[i]);
      vz[i] -= dt2 * odnf * vbar * (ratcom_pz[i] * invMu - wui * vz[i]);
   }
   // iL4 dt/2
   propagate_velocity(dt2, gx, gy, gz);


   // volume iL3
   const double eterm2 = std::exp(vbar * dt);
   lvec1 *= eterm2;
   lvec2 *= eterm2;
   lvec3 *= eterm2;
   set_default_recip_box();


   // position
   // iL2 dt/2
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,mfrac)
   for (int i = 0; i < n; ++i) {
      auto wui = mfrac[i];
      auto coef = std::exp(vbar * wui * dt2);
      xpos[i] *= coef;
      ypos[i] *= coef;
      zpos[i] *= coef;
   }
   // update Ru
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,ratcom_x,ratcom_y,ratcom_z,\
                  mfrac,imol,kmol)
   for (int im = 0; im < nmol; ++im) {
      const int start = imol[im][0];
      const int end = imol[im][1];
      pos_prec rx = 0, ry = 0, rz = 0;
      #pragma acc loop seq
      for (int i = start; i < end; ++i) {
         const int k = kmol[i];
         const auto wui = mfrac[k];
         rx += wui * xpos[k];
         ry += wui * ypos[k];
         rz += wui * zpos[k];
      }
      ratcom_x[im] = rx;
      ratcom_y[im] = ry;
      ratcom_z[im] = rz;
   }
   // iL1 dt, then iL2 dt/2
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,vx,vy,vz,\
                  ratcom_x,ratcom_y,ratcom_z,mfrac,molec)
   for (int i = 0; i < n; ++i) {
      const int im = molec[i];
      const auto wui = mfrac[i];
      const auto coef = std::exp(vbar * wui * dt2);
      xpos[i] = (vx[i] + vbar * (ratcom_x[im] - wui * xpos[i])) * dt * coef;
      ypos[i] = (vy[i] + vbar * (ratcom_y[im] - wui * ypos[i])) * dt * coef;
      zpos[i] = (vz[i] + vbar * (ratcom_z[im] - wui * zpos[i])) * dt * coef;
   }


   // rattle


   energy(vers1);


   // velocity 2/2


   // rattle2


   // thermostat 2/2
   vvlp_hc_pt(dt, R, 2);
}
}
