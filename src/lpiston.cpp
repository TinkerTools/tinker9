#include "lpiston.h"
#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include "nose.h"
#include "random.h"
#include "rattle.h"
#include "tinker_rt.h"
#include "tool/error.h"
#include <cmath>
#include <tinker/detail/bath.hh>
#include <tinker/detail/freeze.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {
double lp_alpha;
double lp_rats1;
double lp_eksum;
double lp_ekin[3][3];
double lp_vir[9];
virial_buffer lp_vir_buf;


//====================================================================//


void lp_atom_kinetic()
{
   kinetic_energy(lp_eksum, lp_ekin, n, mass, vx, vy, vz);
}


void lp_mol_kinetic()
{
   auto& m = rattle_dmol;
   kinetic_energy(lp_eksum, lp_ekin, m.nmol, m.molmass, ratcom_vx, ratcom_vy,
                  ratcom_vz);
}


extern void lp_mol_virial_acc();
extern void lp_mol_virial_cu();
void lp_virial(bool molP)
{
   if (molP) {
      for (int iv = 0; iv < 9; ++iv)
         lp_vir[iv] = 0;
#if TINKER_CUDART
      if (pltfm_config & CU_PLTFM)
         lp_mol_virial_cu();
      else
#endif
         lp_mol_virial_acc();
   } else {
      for (int iv = 0; iv < 9; ++iv)
         lp_vir[iv] = vir[iv];
   }
}


//====================================================================//


double sinh_id(double x)
{
   double y = std::fabs(x);
   if (y <= 1.0e-8)
      return 1.0;
   else
      return std::sinh(y) / y;
}


extern void propagate_pos_raxbv_acc(pos_prec*, pos_prec*, pos_prec*, pos_prec,
                                    pos_prec*, pos_prec*, pos_prec*, pos_prec,
                                    pos_prec*, pos_prec*, pos_prec*);
void propagate_pos_raxbv(pos_prec* r1, pos_prec* r2, pos_prec* r3, pos_prec a,
                         pos_prec* x1, pos_prec* x2, pos_prec* x3, pos_prec b,
                         pos_prec* y1, pos_prec* y2, pos_prec* y3)
{
   propagate_pos_raxbv_acc(r1, r2, r3, a, x1, x2, x3, b, y1, y2, y3);
}


extern void lp_propagate_mol_vel_acc(vel_prec);
void lp_propagate_mol_vel(vel_prec scal)
{
   lp_propagate_mol_vel_acc(scal);
}


void lp_center_of_mass_acc(const pos_prec*, const pos_prec*, const pos_prec*,
                           pos_prec*, pos_prec*, pos_prec*);
void lp_center_of_mass(const pos_prec* atomx, const pos_prec* atomy,
                       const pos_prec* atomz, pos_prec* molx, pos_prec* moly,
                       pos_prec* molz)
{
   static_assert(std::is_same<pos_prec, vel_prec>::value,
                 "pos_prec and vel_prec must be the same type.");
   lp_center_of_mass_acc(atomx, atomy, atomz, molx, moly, molz);
}


//====================================================================//


void lprat_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void lprat_settle_acc(time_prec, const pos_prec*, const pos_prec*,
                      const pos_prec*);
void lprat_ch_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void lprat_methyl_cu(time_prec, const pos_prec*, const pos_prec*,
                     const pos_prec*);
void lprat(time_prec dt, const pos_prec* xold, const pos_prec* yold,
           const pos_prec* zold)
{
   lprat_settle_acc(dt, xold, yold, zold);
   lprat_ch_acc(dt, xold, yold, zold);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      lprat_methyl_cu(dt, xold, yold, zold);
#endif
   lprat_acc(dt, xold, yold, zold);
}


//====================================================================//


static int nrespa, nbaro;
static double rnd;
static double g0, g1, D;
static bool atomT, molT, atomP, molP, constrain, aniso;
static enum { KW_NULL = 0, KW_ATOM, KW_MOL } kw_p;


void vv_lpiston_init()
{
   auto o = stdout;

   // RESPA keyword:    "RESPA-INNER  0.25"
   // Barostat keyword: "BARO-OUTER    6.0"
   const time_prec eps = 1.0 / 1048576; // 2**-20
   nrespa = (int)(time_step / (mdstuf::arespa + eps)) + 1;
   get_kv("VOLUME-TRIAL", nbaro, 1);
   nbaro = std::max(1, nbaro);

   rnd = 0.0;

   // isotropic NPT
   aniso = bath::anisotrop;

   constrain = use_rattle();
   D = 3.0;

   // molecular pressure keyword: "VOLUME-SCALE  MOLECULAR"
   // atomic pressure keyword:    "VOLUME-SCALE     ATOMIC"
   std::string volscale;
   atomP = false, molP = false;
   atomT = false, molT = false;
   if (constrain) {
      kw_p = KW_MOL;
      volscale = "MOLECULAR";
      int val = std::max(rattle_dmol.nmol, 2);
      g1 = D * (val - 1);
      g0 = D * (rattle_dmol.nmol - 1);
      molP = true;
      molT = true;
   } else {
      kw_p = KW_ATOM;
      volscale = "ATOMIC";
      int val = std::max(n, 2);
      g1 = D * (val - 1);
      g0 = D * (n - 1);
      atomP = true;
      atomT = true;
   }
   lp_alpha = 1.0 + D / g1;

   // Nose-Hoover Chain
   double ekt = units::gasconst * bath::kelvin;
   vbar = 0;
   gbar = 0;
   qbar = (g1 + 1) * ekt * bath::taupres * bath::taupres;
   for (int i = 0; i < maxnose; ++i) {
      vnh[i] = 0;
      gnh[i] = 0;
      qnh[i] = ekt * bath::tautemp * bath::tautemp;
   }
   qnh[0] = g0 * ekt * bath::tautemp * bath::tautemp;
   if (nrespa > 1) {
      darray::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);

      // save fast gradients to gx1 etc.
      energy(calc::grad, RESPA_FAST, respa_tsconfig());
      darray::copy(g::q0, n, gx1, gx);
      darray::copy(g::q0, n, gy1, gy);
      darray::copy(g::q0, n, gz1, gz);

      // save slow gradients to gx2 etc.
      energy(calc::grad, RESPA_SLOW, respa_tsconfig());
      darray::copy(g::q0, n, gx2, gx);
      darray::copy(g::q0, n, gy2, gy);
      darray::copy(g::q0, n, gz2, gz);
   }
   // calculate total gradient and atomic virial
   energy(calc::grad + calc::virial);

   if (constrain) {
      darray::allocate(buffer_size(), &lp_vir_buf);
   }

   print(o, "\n");
   print(o, " Friction                     %12.4lf /ps\n", stodyn::friction);
   print(o, " Time constant for const-T    %12.4lf ps\n", bath::tautemp);
   print(o, " Time constant for const-P    %12.4lf ps\n", bath::taupres);
   print(o, " Pressure estimator           %12s\n", volscale);
   print(o, " LP-G                         %12.0lf\n", g0);
   print(o, " LP-G1                        %12.0lf\n", g1);
   print(o, " NRESPA                       %12d\n", nrespa);
   print(o, " NBARO                        %12d\n", nbaro);
   print(o, "\n");
}


void vv_lpiston_destory()
{
   if (constrain)
      darray::deallocate(lp_vir_buf);
   if (nrespa > 1)
      darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);
}


//====================================================================//


static void iso_tp(time_prec dt)
{
   constexpr int ns = 2;
   const double kbt = units::gasconst * bath::kelvin;
   time_prec t, t2, t4, t8, xt4;
   t = dt / ns, t2 = t / 2, t4 = t / 4, t8 = t / 8, xt4 = nbaro * t4;


   double opgxt4, omgxt4, sd;
   opgxt4 = 1.0 + stodyn::friction * xt4;
   omgxt4 = 1.0 - stodyn::friction * xt4;
   sd = std::sqrt(nbaro * dt * 2.0 * stodyn::friction * kbt / qbar) / (4 * ns);


   const virial_prec tr_vir = lp_vir[0] + lp_vir[4] + lp_vir[8];
   const double vol0 = volbox();
   if (atomT) {
      lp_atom_kinetic();
   } else if (molT) {
      lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
      lp_mol_kinetic();
   }


   double eksum0 = lp_eksum, eksum1;
   double velsc0 = 1.0, velsc1;


   double DelP;
   for (int k = 0; k < ns; ++k) {
      eksum1 = eksum0;
      velsc1 = velsc0;

      for (int i = maxnose - 1; i > -1; --i) {
         if (i == 0)
            gnh[i] = (2 * eksum1 - g0 * kbt) / qnh[i];
         else
            gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];

         if (i == maxnose - 1)
            vnh[i] += gnh[i] * t4;
         else {
            double exptm = std::exp(-vnh[i + 1] * t8);
            vnh[i] = (vnh[i] * exptm + gnh[i] * t4) * exptm;
         }
      }

      if (atomP or molP) {
         DelP = lp_alpha * 2 * eksum1 - tr_vir;
         DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
         gbar = DelP / qbar;
         vbar = gbar * xt4 + omgxt4 * vbar + sd * rnd;
      }

      double scal;
      if (atomP or molP)
         scal = std::exp(-t2 * (vnh[0] + lp_alpha * vbar * nbaro));
      else
         scal = std::exp(-t2 * vnh[0]);
      velsc1 *= scal;
      eksum1 *= (scal * scal);

      if (atomP or molP) {
         DelP = lp_alpha * 2 * eksum1 - tr_vir;
         DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
         gbar = DelP / qbar;
         vbar = (gbar * xt4 + vbar + sd * rnd) / opgxt4;
      }

      for (int i = 0; i < maxnose; ++i) {
         if (i == 0)
            gnh[i] = (2 * eksum1 - g0 * kbt) / qnh[i];
         else
            gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];


         if (i == maxnose - 1)
            vnh[i] += gnh[i] * t4;
         else {
            double exptm = std::exp(-vnh[i + 1] * t8);
            vnh[i] = (vnh[i] * exptm + gnh[i] * t4) * exptm;
         }
      }

      eksum0 = eksum1;
      velsc0 = velsc1;
   }


   const double velsc2 = velsc0 * velsc0;
   lp_eksum *= velsc2;
   for (int ii = 0; ii < 3; ++ii)
      for (int jj = 0; jj < 3; ++jj)
         lp_ekin[ii][jj] *= velsc2;
   if (atomT) {
      darray::scale(g::q0, n, velsc0, vx);
      darray::scale(g::q0, n, velsc0, vy);
      darray::scale(g::q0, n, velsc0, vz);
   } else if (molT) {
      lp_propagate_mol_vel(velsc0 - 1.0);
   }
}


void vv_lpiston_npt(int istep, time_prec dt)
{
   bool mid = (nbaro == 1) or ((istep % nbaro) == (nbaro + 1) / 2);
   atomP = false, molP = false;
   if (mid) {
      if (kw_p == KW_ATOM)
         atomP = true;
      else if (kw_p == KW_MOL)
         molP = true;
   }


   int vers1 = rc_flag & calc::vmask;
   if ((istep % inform::iwrite) != 0)
      vers1 &= ~calc::energy;
   if (not mid)
      vers1 &= ~calc::virial;


   time_prec xdt = nbaro * dt, dt2 = 0.5 * dt, xdt2 = 0.5 * xdt;
   time_prec dti = dt / nrespa, dti2 = dt2 / nrespa;
   time_prec xdti = xdt / nrespa, xdti2 = xdt2 / nrespa;


   iso_tp(dt);
   if (nrespa > 1) {
      // gx1: fast; gx2: slow
      propagate_velocity2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
   } else {
      propagate_velocity(dt2, gx, gy, gz);
   }


   if (constrain) {
      darray::copy(g::q0, n, rattle_xold, xpos);
      darray::copy(g::q0, n, rattle_yold, ypos);
      darray::copy(g::q0, n, rattle_zold, zpos);
   }


   virial_prec vir_fast[9] = {0};
   energy_prec esum_f;


   double s = 1.0;
   if (mid) {
      lvec1 *= std::exp(vbar * xdt);
      lvec2 *= std::exp(vbar * xdt);
      lvec3 *= std::exp(vbar * xdt);
      set_default_recip_box();
      double vt2 = vbar * xdti2;
      // double s1b = 1.0;
      // double s1a = std::exp(vt2);
      double s1 = std::exp(vt2) * sinh_id(vt2);
      s = s1;
   }
   for (int ir = 0; ir < nrespa; ++ir) {
      if (mid and kw_p == KW_ATOM) {
         propagate_pos_axbv(std::exp(vbar * xdti), s * dti);
      } else if (mid and kw_p == KW_MOL) {
         lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
         propagate_pos_raxbv(xpos, ypos, zpos, std::exp(vbar * xdti) - 1.0,
                             ratcom_x, ratcom_y, ratcom_z, s * dti, vx, vy, vz);
      } else {
         propagate_pos(dti);
      }

      if (ir < nrespa - 1) {
         copy_pos_to_xyz(false);
         energy(vers1, RESPA_FAST, respa_tsconfig());
         if (vers1 & calc::virial) {
            lp_virial(molP);
            for (int i = 0; i < 9; ++i)
               vir_fast[i] += lp_vir[i];
         }
         propagate_velocity(dti, gx, gy, gz);
      } else {
         if (constrain) {
            lp_rats1 = 1.0 / s;
            lprat(dt, rattle_xold, rattle_yold, rattle_zold);
         }
         copy_pos_to_xyz(true);
      }
   }


   if (nrespa > 1) {
      // fast force
      energy(vers1, RESPA_FAST, respa_tsconfig());
      darray::copy(g::q0, n, gx1, gx);
      darray::copy(g::q0, n, gy1, gy);
      darray::copy(g::q0, n, gz1, gz);
      copy_energy(vers1, &esum_f);
      if (vers1 & calc::virial) {
         lp_virial(molP);
         for (int i = 0; i < 9; ++i)
            vir_fast[i] += lp_vir[i];
      }

      // slow force
      energy(vers1, RESPA_SLOW, respa_tsconfig());
      darray::copy(g::q0, n, gx2, gx);
      darray::copy(g::q0, n, gy2, gy);
      darray::copy(g::q0, n, gz2, gz);
      if (vers1 & calc::energy)
         esum += esum_f;
      if (vers1 & calc::virial) {
         lp_virial(molP);
         for (int iv = 0; iv < 9; ++iv)
            lp_vir[iv] += vir_fast[iv] / nrespa;
      }

      // gx1: fast; gx2: slow
      propagate_velocity2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
   } else {
      energy(vers1);
      if (vers1 & calc::virial)
         lp_virial(molP);
      propagate_velocity(dt2, gx, gy, gz);
   }
   if (mid)
      rnd = normal<double>();
   iso_tp(dt);
   if (constrain)
      rattle2(dt, false);


   if ((istep % inform::iwrite) == 0) {
      lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
      lp_mol_kinetic();
      printf("\n"
             " Current MolKinetic     %12.4lf Kcal/mole of Frame %8d\n",
             lp_eksum, istep / inform::iwrite);
   }
}
}
