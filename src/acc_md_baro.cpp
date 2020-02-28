#define TINKER_ENABLE_LOG 0
#include "log.h"


#include "box.h"
#include "energy.h"
#include "io_fort_str.h"
#include "md.h"
#include "molecule.h"
#include "nblist.h"
#include "random.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/units.hh>


TINKER_NAMESPACE_BEGIN
void monte_carlo_barostat_update_nb_acc(real epot)
{
   fstr_view volscale = bath::volscale;
   double third = 1.0 / 3.0;
   double volmove = bath::volmove;
   double kt = units::gasconst * bath::kelvin;
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
   device_array::copy(PROCEED_NEW_Q, n, x_pmonte, x);
   device_array::copy(PROCEED_NEW_Q, n, y_pmonte, y);
   device_array::copy(PROCEED_NEW_Q, n, z_pmonte, z);


   if (isotropic) {
      double step_rdm = 2 * random<double>() - 1;
      double step = volmove * step_rdm;
      volnew = volold + step;
      double scale = std::pow(volnew / volold, third);
      TINKER_LOG(
         "MC Barostat Volume: random = {:.6f} dV = {:.6f} scale = {:6f}",
         step_rdm, step, scale);


      lvec1 *= scale;
      lvec2 *= scale;
      lvec3 *= scale;
      set_recip_box(lvec1, lvec2, lvec3, recipa, recipb, recipc);


      if (volscale == "MOLECULAR") {
         int nmol = molecule.nmol;
         const auto* imol = molecule.imol;
         const auto* kmol = molecule.kmol;
         const auto* molmass = molecule.molmass;
         real real_scale = scale - 1;
         #pragma acc parallel loop independent async\
                     deviceptr(imol,kmol,mass,molmass,x,y,z)
         for (int i = 0; i < nmol; ++i) {
            real xcm = 0, ycm = 0, zcm = 0;
            int start = imol[i][0];
            int stop = imol[i][1];
            #pragma acc loop seq
            for (int j = start; j < stop; ++j) {
               int k = kmol[j];
               real weigh = mass[k];
               xcm += x[k] * weigh;
               ycm += y[k] * weigh;
               zcm += z[k] * weigh;
            }
            real term = real_scale / molmass[i];
            real xmove = term * xcm;
            real ymove = term * ycm;
            real zmove = term * zcm;
            #pragma acc loop seq
            for (int j = start; j < stop; ++j) {
               int k = kmol[j];
               x[k] += xmove;
               y[k] += ymove;
               z[k] += zmove;
            }
         }
      }
   }

   // get the potential energy and PV work changes for trial move
   nblist_data(rc_evolve);
   energy(calc::energy);
   real enew;
   copy_energy(calc::energy, &enew, nullptr, nullptr, nullptr, nullptr);
   TINKER_LOG("MC Barostat Enew = {:.6f} Eold = {:.6f}", enew, eold);
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
   TINKER_LOG("MC Barostat (kT): dU = {:.6f} dPV = {:.6f} dK = {:.6f}", dpot,
              dpv, dkin);
   TINKER_LOG("MC Barostat Accep. Ratio: {:.6f}; random: {:.6f}; "
              "reject this move if ramdom .gt. Accep. Ratio",
              expterm, exp_rdm);
   if (exp_rdm > expterm) {
      TINKER_LOG("MC Barostat Move Rejected");
      esum = eold;
      set_default_box(boxold);
      device_array::copy(PROCEED_NEW_Q, n, x, x_pmonte);
      device_array::copy(PROCEED_NEW_Q, n, y, y_pmonte);
      device_array::copy(PROCEED_NEW_Q, n, z, z_pmonte);
      nblist_data(rc_evolve);
   } else {
      TINKER_LOG("MC Barostat Move Accepted");
   }
}
TINKER_NAMESPACE_END
