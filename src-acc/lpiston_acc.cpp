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
#include "tinker_rt.h"
#include "tool/error.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/freeze.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {
namespace {
int irespa, ibaro;
double fric, rnd;
double g0, g1, D;
bool atomT, molT, atomP, molP, constrain, aniso;


enum
{
   KW_NULL = 0,
   KW_ATOM,
   KW_MOL
};
int kw_t, kw_p;
}


void vv_lpiston_init()
{
   auto o = stdout;

   // RESPA keywords:    "RESPA-INNER  0.25"
   // Barostat keywords: "BARO-OUTER   6.0"
   const time_prec eps = 1.0 / 1048576; // 2**-20
   irespa = (int)(time_step / (mdstuf::arespa + eps)) + 1;
   double abaro;
   get_kv("BARO-OUTER", abaro, 1.0);
   ibaro = (int)(abaro * 0.001 / (time_step + eps)) + 1;
   if (ibaro % 2 == 0)
      ibaro = ibaro - 1;
   ibaro = std::max(1, ibaro);

   fric = stodyn::friction;
   rnd = 0.0;

   // isotropic NPT
   aniso = bath::anisotrop;

   constrain = use_rattle();
   D = 3.0;

   // molecular pressure keyword: "VOLUME-SCALE  MOLECULAR"
   // atomic pressure keyword:    "VOLUME-SCALE  ATOMIC"
   // default pressure
   fstr_view volscale_f = bath::volscale;
   std::string volscale = volscale_f.trim();
   if (volscale == "ATOMIC") {
      kw_p = KW_ATOM;
   } else if (volscale == "MOLECULAR") {
      kw_p = KW_MOL;
   } else if (constrain) {
      kw_p = KW_MOL;
   } else {
      kw_p = KW_ATOM;
   }
   if (constrain and kw_p == KW_ATOM) {
      TINKER_THROW("NPT RATTLE cannot use atomic volume scaling."
                   " Set \"VOLUME-SCALE  MOLECULAR\" in the key file.");
   }
   if (not constrain and kw_p == KW_MOL) {
      print(o, " No constraints hence setting VOLUME-SCALE to ATOMIC\n");
      kw_p = KW_ATOM;
      volscale = "ATOMIC";
   }
   atomP = false, molP = false;
   if (kw_p == KW_ATOM) {
      volscale = "ATOMIC";
      atomP = true;
      int val = std::max(n, 2);
      g1 = D * (val - 1);
   } else if (kw_p == KW_MOL) {
      volscale = "MOLECULAR";
      molP = true;
      int val = std::max(rattle_dmol.nmol, 2);
      g1 = D * (val - 1);
   }

   // molecular temperature keyword: "VELOCITY-SCALE  MOLECULAR"
   // atomic temperature keyword:    "VELOCITY-SCALE  ATOMIC"
   // default temperature
   std::string velscale_default, velscale;
   if (constrain) {
      velscale_default = "MOLECULAR";
      kw_t = KW_MOL;
   } else {
      velscale_default = "ATOMIC";
      kw_t = KW_ATOM;
   }
   get_kv("VELOCITY-SCALE", velscale, velscale_default);
   if (velscale == "ATOMIC") {
      kw_t = KW_ATOM;
   } else if (velscale == "MOLECULAR") {
      kw_t = KW_MOL;
   }
   atomT = false, molT = false;
   if (kw_t == KW_ATOM) {
      velscale = "ATOMIC";
      atomT = true;
      g0 = D * (n - 1);
      if (constrain) {
         g0 -= freeze::nrat;
      }
   } else if (kw_t == KW_MOL) {
      velscale = "MOLECULAR";
      molT = true;
      g0 = D * (rattle_dmol.nmol - 1);
   }

   lp_alpha = 1.0 + D / g1;


   // Nose-Hoover Chain
   double ekt = units::gasconst * bath::kelvin;
   vbar = 0;
   gbar = 0;
   qbar = (g0 + D) * ekt * bath::taupres * bath::taupres;
   for (int i = 0; i < maxnose; ++i) {
      vnh[i] = 0;
      gnh[i] = 0;
      qnh[i] = ekt * bath::tautemp * bath::tautemp;
   }
   qnh[0] = g0 * ekt * bath::tautemp * bath::tautemp;
   energy(calc::v6);

   if (use_rattle()) {
      darray::allocate(buffer_size(), &lp_vir_buf);
   }

   print(o, "\n");
   print(o, " Friction                        %12.4lf /ps\n", stodyn::friction);
   print(o, " Time constant for the const-T   %12.4lf ps\n", bath::tautemp);
   print(o, " Time constant for the const-P   %12.4lf ps\n", bath::taupres);
   print(o, " Temperature estimator           %12s\n", velscale);
   print(o, " Pressure estimator              %12s\n", volscale);
   print(o, " LP-G                            %12.0lf\n", g0);
   print(o, " LP-G1                           %12.0lf\n", g1);
   print(o, " IRESPA                          %12d\n", irespa);
   print(o, " IBARO                           %12d\n", ibaro);
   print(o, "\n");
}


void vv_lpiston_destory()
{
   if (use_rattle())
      darray::deallocate(lp_vir_buf);
}


void lp_mol_virial_acc() {}


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


void vv_lpiston_npt_acc(int istep, time_prec dt) {}
}
