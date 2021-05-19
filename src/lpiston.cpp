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
double lp_rats2;
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


void lp_atom_virial()
{
   for (int iv = 0; iv < 9; ++iv)
      lp_vir[iv] = vir[iv];
}


extern void lp_mol_virial_acc();
extern void lp_mol_virial_cu();
void lp_mol_virial()
{
   for (int iv = 0; iv < 9; ++iv)
      lp_vir[iv] = 0;
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      lp_mol_virial_cu();
   else
#endif
      lp_mol_virial_acc();
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


void lp_center_of_mass_acc(const pos_prec*, const pos_prec*, const pos_prec*,
                           pos_prec*, pos_prec*, pos_prec*);
void lp_center_of_mass(const pos_prec* atomx, const pos_prec* atomy,
                       const pos_prec* atomz, pos_prec* molx, pos_prec* moly,
                       pos_prec* molz)
{
   static_assert(std::is_same<pos_prec, vel_prec>::value,
                 "pos_prec and vel_pres must be the same type.");
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


void lprat2_acc(time_prec);
void lprat2_settle_acc(time_prec);
void lprat2_ch_acc(time_prec);
void lprat2_methyl_cu(time_prec);
void lprat2(time_prec dt)
{
   lprat2_settle_acc(dt);
   lprat2_ch_acc(dt);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      lprat2_methyl_cu(dt);
#endif
   lprat2_acc(dt);
}


//====================================================================//


static int nrespa, nbaro;
static double fric, rnd;
static double g0, g1, D;
static bool atomT, molT, atomP, molP, constrain, aniso;
static enum { KW_NULL = 0, KW_ATOM, KW_MOL } kw_t, kw_p;


void vv_lpiston_init()
{
   auto o = stdout;

   // RESPA keywords:    "RESPA-INNER  0.25"
   // Barostat keywords: "BARO-OUTER   6.0"
   const time_prec eps = 1.0 / 1048576; // 2**-20
   nrespa = (int)(time_step / (mdstuf::arespa + eps)) + 1;
   double abaro;
   get_kv("BARO-OUTER", abaro, 1.0);
   nbaro = (int)(abaro * 0.001 / (time_step + eps)) + 1;
   if (nbaro % 2 == 0)
      nbaro = nbaro - 1;
   nbaro = std::max(1, nbaro);

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
   print(o, " NRESPA                          %12d\n", nrespa);
   print(o, " NBARO                           %12d\n", nbaro);
   print(o, "\n");
}


void vv_lpiston_destory()
{
   if (use_rattle())
      darray::deallocate(lp_vir_buf);
   if (nrespa > 1)
      darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);
}


//====================================================================//


static void iso_tp(int half, time_prec dt)
{
   constexpr int ns = 2;
   const double kbt = units::gasconst * bath::kelvin;
   time_prec DT = nbaro * dt, dt2 = 0.5 * dt, DT2 = 0.5 * DT;
   time_prec h, h2, h4, H, H2;
   h = dt2 / ns, h2 = 0.5 * h, h4 = 0.25 * h;
   H = DT2 / ns, H2 = 0.5 * H;
   double opgh2, omgh2, sd, alpha;
   opgh2 = 1.0 + fric * H2;
   omgh2 = 1.0 - fric * H2;
   sd = sqrt(2.0 * kbt * fric * DT / qbar) / (4 * n);
   alpha = lp_alpha;
}


void vv_lpiston_npt(int istep, time_prec dt)
{
   time_prec DT = nbaro * dt, dt2 = 0.5 * dt, DT2 = 0.5 * DT;
   time_prec dti = dt / nrespa, dti2 = dt2 / nrespa;
   time_prec DTi = DT / nrespa, DTi2 = DT2 / nrespa;


   bool mid = (istep % nbaro) == (nbaro + 1) / 2;
   lp_rats1 = 1.0, lp_rats2 = 1.0;


   iso_tp(1, dt);
   if (nrespa > 1) {
      // gx1: fast; gx2: slow
      propagate_velocity2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
   } else {
      propagate_velocity(dt, gx, gy, gz);
   }


   if (constrain) {
      darray::copy(g::q0, n, rattle_xold, xpos);
      darray::copy(g::q0, n, rattle_yold, ypos);
      darray::copy(g::q0, n, rattle_zold, zpos);
   }


   double s;
   if (mid and (atomP or molP)) {
      lvec1 *= std::exp(vbar * DT);
      lvec2 *= std::exp(vbar * DT);
      lvec3 *= std::exp(vbar * DT);
      set_default_recip_box();
      double vt2 = vbar * DTi2;
      // double s1b = 1.0;
      // double s1a = std::exp(vt2);
      double s1 = std::exp(vt2) * sinh_id(vt2);
      s = s1;
      lp_rats1 = 1.0 / s;
   }
   for (int ir = 0; ir < nrespa; ++ir) {
      if (mid and atomP) {
         propagate_pos_axbv(std::exp(vbar * DTi), s * dti);
      } else if (mid and molP) {
         lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
         #pragma acc parallel loop independent async\
                 deviceptr(vx,vy,vz,xpos,ypos,zpos,ratcom_x,ratcom_y,ratcom_z)
         for (int i = 0; i < n; ++i) {
            double c1 = s * dti, c2 = std::exp(vbar * DTi) - 1.0;
            xpos[i] += c1 * vx[i] + c2 * ratcom_x[i];
            ypos[i] += c1 * vy[i] + c2 * ratcom_y[i];
            zpos[i] += c1 * vz[i] + c2 * ratcom_z[i];
         }
      } else {
         propagate_pos(dti);
      }

      if (ir < nrespa - 1) {
         energy(calc::grad, RESPA_FAST, respa_tsconfig());
         propagate_velocity(dti, gx, gy, gz);
      }
   }
   if (constrain) {
      assert(not atomP);
      lprat(dt, rattle_xold, rattle_yold, rattle_zold);
   }


   rnd = normal<double>();
   atomP = false, molP = false;
   if (istep % nbaro == 0) {
      if (kw_p == KW_ATOM)
         atomP = true;
      if (kw_p == KW_MOL)
         molP = true;
   }


   int vers1 = rc_flag & calc::vmask;
   if ((istep % inform::iwrite) != 0)
      vers1 &= ~calc::energy;
   if ((istep % nbaro) != 0)
      vers1 &= ~calc::virial;


   if (nrespa > 1) {
      // fast force
      energy(vers1, RESPA_FAST, respa_tsconfig());
      darray::copy(g::q0, n, gx1, gx);
      darray::copy(g::q0, n, gy1, gy);
      darray::copy(g::q0, n, gz1, gz);
      energy_prec esum_f;
      virial_prec vir_f[9];
      copy_energy(vers1, &esum_f);
      copy_virial(vers1, vir_f);


      // slow force
      energy(vers1, RESPA_SLOW, respa_tsconfig());
      darray::copy(g::q0, n, gx2, gx);
      darray::copy(g::q0, n, gy2, gy);
      darray::copy(g::q0, n, gz2, gz);
      if (vers1 & calc::energy)
         esum += esum_f;
      if (vers1 & calc::virial) {
         for (int iv = 0; iv < 9; ++iv)
            vir[iv] += vir_f[iv];
      }


      // gx1: fast; gx2: slow
      propagate_velocity2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
   } else {
      energy(vers1);
      propagate_velocity(dt, gx, gy, gz);
   }
   iso_tp(2, dt);
   if (constrain)
      lprat2(dt);
}
}
