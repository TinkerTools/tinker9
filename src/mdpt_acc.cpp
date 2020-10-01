#define TINKER_ENABLE_LOG 0
#include "tool/log.h"


#include "box.h"
#include "energy.h"
#include "glob.molecule.h"
#include "mathfunc.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "mdpt.h"
#include "nblist.h"
#include "random.h"
#include "tool/io_fort_str.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/bound.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>


namespace tinker {
void kinetic_explicit_acc(T_prec& temp_out, energy_prec& eksum_out,
                          energy_prec (&ekin_out)[3][3], const vel_prec* vx,
                          const vel_prec* vy, const vel_prec* vz)
{
   const energy_prec ekcal_inv = 1.0 / units::ekcal;
   energy_prec exx = 0;
   energy_prec eyy = 0;
   energy_prec ezz = 0;
   energy_prec exy = 0;
   energy_prec eyz = 0;
   energy_prec ezx = 0;
   #pragma acc parallel loop independent async\
               deviceptr(mass,vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      energy_prec term = 0.5f * mass[i] * ekcal_inv;
      exx += term * vx[i] * vx[i];
      eyy += term * vy[i] * vy[i];
      ezz += term * vz[i] * vz[i];
      exy += term * vx[i] * vy[i];
      eyz += term * vy[i] * vz[i];
      ezx += term * vz[i] * vx[i];
   }
   energy_prec eksum_local = exx + eyy + ezz;
   T_prec temp_local = 2 * eksum_local / (mdstuf::nfree * units::gasconst);


   ekin_out[0][0] = exx;
   ekin_out[0][1] = exy;
   ekin_out[0][2] = ezx;
   ekin_out[1][0] = exy;
   ekin_out[1][1] = eyy;
   ekin_out[1][2] = eyz;
   ekin_out[2][0] = ezx;
   ekin_out[2][1] = eyz;
   ekin_out[2][2] = ezz;
   eksum_out = eksum_local;
   temp_out = temp_local;
}


//====================================================================//


void bussi_thermostat_acc(time_prec dt_prec, T_prec temp_prec)
{
   double dt = dt_prec;
   double temp = temp_prec;


   double tautemp = bath::tautemp;
   double kelvin = bath::kelvin;
   int nfree = mdstuf::nfree;
   double& eta = bath::eta;


   if (temp == 0)
      temp = 0.1;


   double c = std::exp(-dt / tautemp);
   double d = (1 - c) * (kelvin / temp) / nfree;
   double r = normal<double>();
   double s = chi_squared<double>(nfree - 1);
   double scale = c + (s + r * r) * d + 2 * r * std::sqrt(c * d);
   scale = std::sqrt(scale);
   if (r + std::sqrt(c / d) < 0)
      scale = -scale;
   eta *= scale;


   vel_prec sc = scale;
   #pragma acc parallel loop independent async deviceptr(vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      vx[i] *= sc;
      vy[i] *= sc;
      vz[i] *= sc;
   }
}


//====================================================================//


void monte_carlo_barostat_acc(energy_prec epot, T_prec temp)
{
   if (not bound::use_bounds)
      return;
   if (thermostat != NONE_THERMOSTAT)
      temp = bath::kelvin;


   fstr_view volscale = bath::volscale;
   double third = 1.0 / 3.0;
   double volmove = bath::volmove;
   double kt = units::gasconst * temp;
   bool isotropic = true;
   // double aniso_rdm = random<double>();
   // if (bath::anisotrop && aniso_rdm > 0.5)
   //    isotropic = false;


   // save the system state prior to trial box size change
   Box boxold;
   get_default_box(boxold);
   double volold = volbox();
   double volnew = 0;
   double eold = epot;
   darray::copy(PROCEED_NEW_Q, n, x_pmonte, xpos);
   darray::copy(PROCEED_NEW_Q, n, y_pmonte, ypos);
   darray::copy(PROCEED_NEW_Q, n, z_pmonte, zpos);
   darray::copy(PROCEED_NEW_Q, n, vx_pmonte, vx);
   darray::copy(PROCEED_NEW_Q, n, vy_pmonte, vy);
   darray::copy(PROCEED_NEW_Q, n, vz_pmonte, vz);


   if (isotropic) {
      double step_rdm = 2 * random<double>() - 1;
      double step = volmove * step_rdm;
      volnew = volold + step;
      double scale = std::pow(volnew / volold, third);
      TINKER_LOG("MC Barostat Isotropic: random = %.6f dV = %.6f scale = %.6f",
                 step_rdm, step, scale);


      lvec1 *= scale;
      lvec2 *= scale;
      lvec3 *= scale;
      set_default_recip_box();


      if (volscale == "MOLECULAR") {
         int nmol = molecule.nmol;
         const auto* imol = molecule.imol;
         const auto* kmol = molecule.kmol;
         const auto* molmass = molecule.molmass;
         pos_prec pos_scale = scale - 1;
         #pragma acc parallel loop independent async\
                     deviceptr(imol,kmol,mass,molmass,xpos,ypos,zpos,vx,vy,vz)
         for (int i = 0; i < nmol; ++i) {
            pos_prec xcm = 0, ycm = 0, zcm = 0;
            vel_prec vxcm = 0, vycm = 0, vzcm = 0;
            int start = imol[i][0];
            int stop = imol[i][1];
            #pragma acc loop seq
            for (int j = start; j < stop; ++j) {
               int k = kmol[j];
               mass_prec weigh = mass[k];
               xcm += xpos[k] * weigh;
               ycm += ypos[k] * weigh;
               zcm += zpos[k] * weigh;
               vxcm += vx[k] * weigh;
               vycm += vy[k] * weigh;
               vzcm += vz[k] * weigh;
            }
            pos_prec term = pos_scale / molmass[i];
            pos_prec xmove, ymove, zmove;
            vel_prec vxmove, vymove, vzmove;
            xmove = term * xcm;
            ymove = term * ycm;
            zmove = term * zcm;
            vxmove = term * vxcm;
            vymove = term * vycm;
            vzmove = term * vzcm;
            #pragma acc loop seq
            for (int j = start; j < stop; ++j) {
               int k = kmol[j];
               xpos[k] += xmove;
               ypos[k] += ymove;
               zpos[k] += zmove;
               vx[k] -= vxmove;
               vy[k] -= vymove;
               vz[k] -= vzmove;
            }
         }
         copy_pos_to_xyz();
      }
   }

   // get the potential energy and PV work changes for trial move
   refresh_neighbors();
   energy(calc::energy);
   energy_prec enew;
   copy_energy(calc::energy, &enew);
   TINKER_LOG("MC Barostat Enew = %.6f Eold = %.6f", enew, eold);
   double dpot = enew - eold;
   double dpv = bath::atmsph * (volnew - volold) / units::prescon;


   // estimate the kinetic energy change as an ideal gas term
   double dkin = 0;
   if (volscale == "MOLECULAR") {
      dkin = molecule.nmol * kt * std::log(volold / volnew);
   }


   // acceptance ratio from Epot change, Ekin change and PV work
   double term = -(dpot + dpv + dkin) / kt;
   double expterm = std::exp(term);


   // reject the step, and restore values prior to trial change
   double exp_rdm = random<double>();
   TINKER_LOG("MC Barostat (kT): dU = %.6f dPV = %.6f dK = %.6f", dpot, dpv,
              dkin);
   TINKER_LOG("MC Barostat Accep. Ratio: %.6f; random: %.6f; "
              "reject this move if ramdom .gt. Accep. Ratio",
              expterm, exp_rdm);
   if (exp_rdm > expterm) {
      TINKER_LOG("MC Barostat Move Rejected");
      esum = eold;
      set_default_box(boxold);
      darray::copy(PROCEED_NEW_Q, n, xpos, x_pmonte);
      darray::copy(PROCEED_NEW_Q, n, ypos, y_pmonte);
      darray::copy(PROCEED_NEW_Q, n, zpos, z_pmonte);
      darray::copy(PROCEED_NEW_Q, n, vx, vx_pmonte);
      darray::copy(PROCEED_NEW_Q, n, vy, vy_pmonte);
      darray::copy(PROCEED_NEW_Q, n, vz, vz_pmonte);
      copy_pos_to_xyz();
      refresh_neighbors();
   } else {
#if TINKER_ENABLE_LOG
      Box p;
      get_default_box(p);
      double xbox, ybox, zbox, a_deg, b_deg, c_deg;
      get_box_axes_angles(p, xbox, ybox, zbox, a_deg, b_deg, c_deg);
      TINKER_LOG("MC Barostat Move Accepted; New box"_s + 6 * "%12.6f"_s, xbox,
                 ybox, zbox, a_deg, b_deg, c_deg);
#endif
   }
}
}
