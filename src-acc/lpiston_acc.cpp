#include "add.h"
#include "box.h"
#include "energy.h"
#include "lpiston.h"
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
void lp_molpressure_acc(double alpha, double& val)
{
   int nmol = rattle_dmol.nmol;
   auto* molmass = rattle_dmol.molmass;
   auto* imol = rattle_dmol.imol;
   auto* kmol = rattle_dmol.kmol;


   double sum = 0;
   #pragma acc parallel loop independent async\
               copy(sum) reduction(+:sum)\
               deviceptr(molmass,imol,kmol,ratcom_x,ratcom_y,ratcom_z,\
                  mass,xpos,ypos,zpos,vx,vy,vz,gx,gy,gz)
   for (int im = 0; im < nmol; ++im) {
      double mvir = 0.0;
      double igx, igy, igz;
      double mgx = 0, mgy = 0, mgz = 0;
      vel_prec px = 0, py = 0, pz = 0;
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


      sum += (alpha * mv2 - mvir);
   }
   #pragma acc wait


   virial_prec atomic_vir = vir[0] + vir[4] + vir[8];
   val = sum + atomic_vir;
}


namespace {
void ratcom_p()
{
   const int nmol = rattle_dmol.nmol;
   const auto* imol = rattle_dmol.imol;
   const auto* kmol = rattle_dmol.kmol;
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,ratcom_x,ratcom_y,ratcom_z,\
                  vx,vy,vz,ratcom_px,ratcom_py,ratcom_pz,\
                  mass,imol,kmol)
   for (int im = 0; im < nmol; ++im) {
      int start = imol[im][0];
      int end = imol[im][1];
      vel_prec ptx = 0, pty = 0, ptz = 0;
      #pragma acc loop seq
      for (int i = start; i < end; ++i) {
         int k = kmol[i];
         auto massk = mass[k];
         ptx += massk * vx[k];
         pty += massk * vy[k];
         ptz += massk * vz[k];
      }
      ratcom_px[im] = ptx;
      ratcom_py[im] = pty;
      ratcom_pz[im] = ptz;
   }
}


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
   const double odnf = lp_alpha;
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


   double vol0 = volbox();
   lp_molpressure(odnf, lp_molpres);
   gbar = lp_molpres - D * vol0 * bath::atmsph / units::prescon;
   gbar /= qbar;
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


void vv_lpiston_hc_acc(int istep, time_prec dt)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(istep % inform::iwrite);
   if (!save)
      vers1 &= ~calc::energy;


   const time_prec dt2 = 0.5 * dt;
   const double odnf = lp_alpha;
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
      auto coef = exp(-dt2 * (vnh0 + odnf * vbar * wui));
      vx[i] *= coef;
      vy[i] *= coef;
      vz[i] *= coef;
   }
   // iL5 dt/2
   ratcom_p();
   #pragma acc parallel loop independent async\
               deviceptr(vx,vy,vz,ratcom_px,ratcom_py,ratcom_pz,\
                  mfrac,molmass,molec)
   for (int i = 0; i < n; ++i) {
      auto im = molec[i];
      auto invMu = 1.0 / molmass[im];
      auto wui = mfrac[i];
      vx[i] -= dt2 * odnf * vbar * (ratcom_px[im] * invMu - wui * vx[i]);
      vy[i] -= dt2 * odnf * vbar * (ratcom_py[im] * invMu - wui * vy[i]);
      vz[i] -= dt2 * odnf * vbar * (ratcom_pz[im] * invMu - wui * vz[i]);
   }
   // iL4 dt/2
   propagate_velocity(dt2, gx, gy, gz);


   // volume iL3
   const double eterm2 = std::exp(vbar * dt);
   lvec1 *= eterm2;
   lvec2 *= eterm2;
   lvec3 *= eterm2;
   set_default_recip_box();


   darray::copy(g::q0, n, rattle_xold, xpos);
   darray::copy(g::q0, n, rattle_yold, ypos);
   darray::copy(g::q0, n, rattle_zold, zpos);


   // position
   // iL2 dt/2, then center-of-mass
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,ratcom_x,ratcom_y,ratcom_z,\
                  mfrac,imol,kmol)
   for (int im = 0; im < nmol; ++im) {
      int start = imol[im][0];
      int end = imol[im][1];
      pos_prec rtx = 0, rty = 0, rtz = 0;
      #pragma acc loop seq
      for (int i = start; i < end; ++i) {
         int k = kmol[i];
         auto wui = mfrac[k];
         auto coef = exp(vbar * wui * dt2);

         pos_prec xi, yi, zi;
         xi = xpos[k] * coef;
         yi = ypos[k] * coef;
         zi = zpos[k] * coef;
         xpos[k] = xi;
         ypos[k] = yi;
         zpos[k] = zi;

         rtx += wui * xi;
         rty += wui * yi;
         rtz += wui * zi;
      }
      ratcom_x[im] = rtx;
      ratcom_y[im] = rty;
      ratcom_z[im] = rtz;
   }


   // iL1 dt, then iL2 dt/2
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,vx,vy,vz,\
                  ratcom_x,ratcom_y,ratcom_z,mfrac,molec)
   for (int i = 0; i < n; ++i) {
      const int im = molec[i];
      const auto wui = mfrac[i];
      const auto coef = exp(vbar * wui * dt2);
      pos_prec xi, yi, zi;
      xi = xpos[i];
      yi = ypos[i];
      zi = zpos[i];
      xi = xi + dt * (vx[i] + vbar * (ratcom_x[im] - wui * xi));
      yi = yi + dt * (vy[i] + vbar * (ratcom_y[im] - wui * yi));
      zi = zi + dt * (vz[i] + vbar * (ratcom_z[im] - wui * zi));
      xpos[i] = xi * coef;
      ypos[i] = yi * coef;
      zpos[i] = zi * coef;
   }


   // rattle
   lprat(dt, rattle_xold, rattle_yold, rattle_zold);
   copy_pos_to_xyz(true);


   // update gradient
   energy(vers1);


   // velocity 2/2
   // iL4 dt/2
   propagate_velocity(dt2, gx, gy, gz);
   // iL5 dt/2
   ratcom_p();
   #pragma acc parallel loop independent async\
               deviceptr(vx,vy,vz,ratcom_px,ratcom_py,ratcom_pz,\
                  mfrac,molmass,molec)
   for (int i = 0; i < n; ++i) {
      auto im = molec[i];
      auto invMu = 1.0 / molmass[im];
      auto wui = mfrac[i];
      vx[i] -= dt2 * odnf * vbar * (ratcom_px[im] * invMu - wui * vx[i]);
      vy[i] -= dt2 * odnf * vbar * (ratcom_py[im] * invMu - wui * vy[i]);
      vz[i] -= dt2 * odnf * vbar * (ratcom_pz[im] * invMu - wui * vz[i]);
   }
   // iL6 dt/2
   #pragma acc parallel loop independent async\
               deviceptr(vx,vy,vz,mfrac)
   for (int i = 0; i < n; ++i) {
      auto wui = mfrac[i];
      auto coef = exp(-dt2 * (vnh0 + odnf * vbar * wui));
      vx[i] *= coef;
      vy[i] *= coef;
      vz[i] *= coef;
   }


   // rattle2
   lprat2(dt);


   // thermostat 2/2
   vvlp_hc_pt(dt, R, 2);
}
}
