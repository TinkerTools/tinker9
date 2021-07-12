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
#include "tool/orthomatrix.h"
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


double vbar_matrix[3][3];
double vbar_eigen[3];
double vbar_ortho[3][3];


//====================================================================//


extern void lp_matvec_acc(int len, char transpose, double mat[3][3],
                          pos_prec* ax, pos_prec* ay, pos_prec* az);
void lp_matvec(int len, char transpose, double mat[3][3], pos_prec* ax,
               pos_prec* ay, pos_prec* az)
{
   lp_matvec_acc(len, transpose, mat, ax, ay, az);
}


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


extern void propagate_pos_raxbv_acc(

   pos_prec*, pos_prec*, pos_prec*,

   const int*, double, pos_prec*, pos_prec*, pos_prec*,

   double, pos_prec*, pos_prec*, pos_prec*);
void propagate_pos_raxbv(double a, double b)
{
   propagate_pos_raxbv_acc(xpos, ypos, zpos,

                           rattle_dmol.molecule,

                           a, ratcom_x, ratcom_y, ratcom_z,

                           b, vx, vy, vz);
}


extern void propagate_pos_raxbv_aniso_acc(

   pos_prec*, pos_prec*, pos_prec*,

   const int*, double[3][3], pos_prec*, pos_prec*, pos_prec*,

   double[3][3], pos_prec*, pos_prec*, pos_prec*);
void propagate_pos_raxbv_aniso(double a[3][3], double b[3][3])
{
   propagate_pos_raxbv_aniso_acc(xpos, ypos, zpos,

                                 rattle_dmol.molecule,

                                 a, ratcom_x, ratcom_y, ratcom_z,

                                 b, vx, vy, vz);
}


extern void propagate_pos_axbv_aniso_acc(double[3][3], double[3][3]);
void propagate_pos_axbv_aniso(double a[3][3], double b[3][3])
{
   propagate_pos_axbv_aniso_acc(a, b);
}


extern void lp_propagate_mol_vel_acc(vel_prec);
void lp_propagate_mol_vel(vel_prec scal)
{
   lp_propagate_mol_vel_acc(scal);
}


extern void lp_propagate_mol_vel_aniso_acc(vel_prec scal[3][3]);
void lp_propagate_mol_vel_aniso(vel_prec scal[3][3])
{
   lp_propagate_mol_vel_aniso_acc(scal);
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
static double rnd, rnd_matrix[3][3];
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
   std::string itg;
   get_kv("INTEGRATOR", itg, "");
   if (itg != "RESPA")
      nrespa = 1;
   get_kv("VOLUME-TRIAL", nbaro, 1);
   nbaro = std::max(1, nbaro);

   get_kv("FRICTION", stodyn::friction, 20.0);

   rnd = 0.0;
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         rnd_matrix[i][j] = 0.0;

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

   // isotropic NPT
   aniso = bath::anisotrop;
   if (aniso) {
      qbar /= D;
      vbar_matrix[0][0] = 0;
      vbar_matrix[0][1] = 0;
      vbar_matrix[0][2] = 0;
      vbar_matrix[1][0] = 0;
      vbar_matrix[1][1] = 0;
      vbar_matrix[1][2] = 0;
      vbar_matrix[2][0] = 0;
      vbar_matrix[2][1] = 0;
      vbar_matrix[2][2] = 0;
      vbar_eigen[0] = 0;
      vbar_eigen[1] = 0;
      vbar_eigen[2] = 0;
      vbar_ortho[0][0] = 1.0;
      vbar_ortho[0][1] = 0;
      vbar_ortho[0][2] = 0;
      vbar_ortho[1][0] = 0;
      vbar_ortho[1][1] = 1.0;
      vbar_ortho[1][2] = 0;
      vbar_ortho[2][0] = 0;
      vbar_ortho[2][1] = 0;
      vbar_ortho[2][2] = 1.0;
   }

   print(o, "\n");
   print(o, " Friction                     %12.4lf /ps\n", stodyn::friction);
   print(o, " Time constant for const-T    %12.4lf ps\n", bath::tautemp);
   print(o, " Time constant for const-P    %12.4lf ps\n", bath::taupres);
   print(o, " Pressure estimator           %12s\n", volscale);
   print(o, " Kinetic DOF                  %12.0lf\n", g0);
   print(o, " Pressure DOF                 %12.0lf\n", g1);
   print(o, " NRESPA                       %12d\n", nrespa);
   print(o, " NBARO                        %12d\n", nbaro);
   if (aniso)
      print(o, " ANISOTROPIC\n");
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
   for (int k = 0; k < ns; ++k) {
      double DelP;

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


static void lpiston_npt_aniso(int, time_prec);
void vv_lpiston_npt(int istep, time_prec dt)
{
   if (aniso) {
      lpiston_npt_aniso(istep, dt);
      return;
   }


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
         propagate_pos_raxbv(std::exp(vbar * xdti) - 1.0, s * dti);
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


   if (molT and (istep % inform::iwrite) == 0) {
      lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
      lp_mol_kinetic();
      printf("\n"
             " Current MolKinetic     %12.4lf Kcal/mole of Frame %8d\n",
             lp_eksum, istep / inform::iwrite);
   }
}


//====================================================================//


template <class T>
static void prtmat(T m[3][3], const char* flag = nullptr)
{
   if (flag)
      printf(" %s\n", flag);
   for (int i = 0; i < 3; ++i) {
      if (i == 0)
         printf(" [[");
      else
         printf("  [");
      printf("%16.4e,%16.4e,%16.4e]", m[i][0], m[i][1], m[i][2]);
      if (i == 0 or i == 1)
         printf(",\n");
      if (i == 2)
         printf("]\n");
   }
   printf("\n");
}

template <class T>
static void prtvec(T m[3], const char* flag = nullptr)
{
   if (flag)
      printf(" %s\n", flag);
   printf("  [%16.4e,%16.4e,%16.4e]\n", m[0], m[1], m[2]);
   printf("\n");
}


static void iso_tp_aniso(time_prec dt)
{
   // constexpr int ns = 2;
   constexpr int ns = 1;
   const double kbt = units::gasconst * bath::kelvin;
   time_prec t, t2, t4, t8, xt4;
   t = dt / ns, t2 = t / 2, t4 = t / 4, t8 = t / 8, xt4 = nbaro * t4;


   double opgxt4, omgxt4, sd;
   opgxt4 = 1.0 + stodyn::friction * xt4;
   omgxt4 = 1.0 - stodyn::friction * xt4;
   sd = std::sqrt(nbaro * dt * 2.0 * stodyn::friction * kbt / qbar) / (4 * ns);


   const double vol0 = volbox();
   if (atomT) {
      lp_atom_kinetic();
   } else if (molT) {
      lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
      lp_mol_kinetic();
   }


   double eksum1 = lp_eksum,
          ekin0[3][3] = {{lp_ekin[0][0], lp_ekin[0][1], lp_ekin[0][2]},
                         {lp_ekin[1][0], lp_ekin[1][1], lp_ekin[1][2]},
                         {lp_ekin[2][0], lp_ekin[2][1], lp_ekin[2][2]}};
   double velsc0[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}},
          velsc1[3][3];
   for (int k = 0; k < ns; ++k) {
      double DelP[3][3];


      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            velsc1[i][j] *= velsc0[i][j];


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
         for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
               DelP[i][j] = 2 * ekin0[i][j] - lp_vir[3 * i + j];
               if (j == i) {
                  DelP[i][i] +=
                     2 * eksum1 / g1 - vol0 * bath::atmsph / units::prescon;
               }
               DelP[i][j] /= qbar;
               vbar_matrix[i][j] = DelP[i][j] * xt4 +
                  omgxt4 * vbar_matrix[i][j] + sd * rnd_matrix[i][j];
            }
            for (int j = 0; j < i; ++j) {
               vbar_matrix[i][j] = vbar_matrix[j][i];
            }
         }
         SymmMatrix::solve(vbar_matrix, vbar_ortho, vbar_eigen);
         prtmat(vbar_matrix, "vbar");
         prtmat(vbar_ortho, "vbar_ortho");
         prtvec(vbar_eigen, "vbar_eigen");
      }


      if (atomP or molP) {
         double tr;
         tr = (vbar_matrix[0][0] + vbar_matrix[1][1] + vbar_matrix[2][2]) / g1;
         double diag[3] = {
            std::exp(-t2 * (vnh[0] + tr * vbar_eigen[0] * nbaro)),
            std::exp(-t2 * (vnh[0] + tr * vbar_eigen[1] * nbaro)),
            std::exp(-t2 * (vnh[0] + tr * vbar_eigen[2] * nbaro))};
         // V(3) = O(3,3) DIAG(3) O^T(3,3) V(3)
         // Scal(3,3) = O(3,3) DIAG(3) O^T(3,3)
         double scal[3][3];
         SymmMatrix::ODOt(scal, vbar_ortho, diag);
         matmul3(velsc1, scal);

         // kinetic tensor = ODO^T[3][3] K[3][3] ODO^T[3][3]
         double a = scal[0][0], b = scal[0][1], c = scal[0][2];
         double d = scal[1][1], e = scal[1][2], f = scal[2][2];
         double u = ekin0[0][0], v = ekin0[0][1], w = ekin0[0][2];
         double x = ekin0[1][1], y = ekin0[1][2], z = ekin0[2][2];

         // [[u,v,w]    [[a,b,c]  [[u,v,w]  [[a,b,c]
         //  [v,x,y]  =  [b,d,e]   [v,x,y]   [b,d,e]
         //  [w,y,z]]    [c,e,f]]  [w,y,z]]  [c,e,f]]

         // clang-format off
         ekin0[0][0] = a*a*u + 2*a*b*v + 2*a*c*w + b*b*x + 2*b*c*y + c*c*z;
         ekin0[1][1] = b*b*u + 2*b*d*v + 2*b*e*w + d*d*x + 2*d*e*y + e*e*z;
         ekin0[2][2] = c*c*u + 2*c*e*v + 2*c*f*w + e*e*x + 2*e*f*y + f*f*z;
         ekin0[0][1] = a*b*u + a*d*v + a*e*w + b*b*v + b*c*w + b*d*x + b*e*y + c*d*y + c*e*z;
         ekin0[0][2] = a*c*u + a*e*v + a*f*w + b*c*v + b*e*x + b*f*y + c*c*w + c*e*y + c*f*z;
         // ekin0[1][0] = a*b*u + a*d*v + a*e*w + b*b*v + b*c*w + b*d*x + b*e*y + c*d*y + c*e*z;
         ekin0[1][0] = ekin0[0][1];
         ekin0[1][2] = b*c*u + b*e*v + b*f*w + c*d*v + c*e*w + d*e*x + d*f*y + e*e*y + e*f*z;
         // ekin0[2][0] = a*c*u + a*e*v + a*f*w + b*c*v + b*e*x + b*f*y + c*c*w + c*e*y + c*f*z;
         // ekin0[2][1] = b*c*u + b*e*v + b*f*w + c*d*v + c*e*w + d*e*x + d*f*y + e*e*y + e*f*z;
         ekin0[2][0] = ekin0[0][2];
         ekin0[2][1] = ekin0[1][2];
         // clang-format on
         eksum1 = ekin0[0][0] + ekin0[1][1] + ekin0[2][2];
      } else {
         double scal = std::exp(-t2 * vnh[0]);
         scal *= scal;
         for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
               ekin0[i][j] *= scal;
         eksum1 *= scal;
      }


      if (atomP or molP) {
         for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
               DelP[i][j] = 2 * ekin0[i][j] - lp_vir[3 * i + j];
               if (j == i) {
                  DelP[i][i] +=
                     2 * eksum1 / g1 - vol0 * bath::atmsph / units::prescon;
               }
               DelP[i][j] /= qbar;
               vbar_matrix[i][j] = (DelP[i][j] * xt4 + vbar_matrix[i][j] +
                                    sd * rnd_matrix[i][j]) /
                  opgxt4;
            }
            for (int j = 0; j < i; ++j) {
               vbar_matrix[i][j] = vbar_matrix[j][i];
            }
         }
         SymmMatrix::solve(vbar_matrix, vbar_ortho, vbar_eigen);
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


      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            velsc0[i][j] *= velsc1[i][j];
   }


   for (int ii = 0; ii < 3; ++ii)
      for (int jj = 0; jj < 3; ++jj)
         lp_ekin[ii][jj] = ekin0[ii][jj];
   lp_eksum = lp_ekin[0][0] + lp_ekin[1][1] + lp_ekin[2][2];
   if (atomT) {
      lp_matvec(n, 'N', velsc0, vx, vy, vz);
   } else if (molT) {
      vel_prec scal[3][3];
      for (int i = 0; i < 3; ++i) {
         for (int j = 0; j < 3; ++j) {
            scal[i][j] = velsc0[i][j];
         }
         scal[i][i] -= 1.0;
      }
      lp_propagate_mol_vel_aniso(scal);
   }
}


static void lpiston_npt_aniso(int istep, time_prec dt)
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


   iso_tp_aniso(dt);
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


   double s[3] = {1.0, 1.0, 1.0};
   if (mid) {
      double diag[3] = {std::exp(vbar_eigen[0] * xdt),
                        std::exp(vbar_eigen[1] * xdt),
                        std::exp(vbar_eigen[2] * xdt)},
             scal[3][3];
      SymmMatrix::ODOt(scal, vbar_ortho, diag);
      double h0[3][3] = {{lvec1.x, lvec1.y, lvec1.z},
                         {lvec2.x, lvec2.y, lvec2.z},
                         {lvec3.x, lvec3.y, lvec3.z}};
      matmul3(h0, scal);
      lvec1.x = h0[0][0], lvec1.y = h0[0][1], lvec1.z = h0[0][2];
      lvec2.x = h0[1][0], lvec2.y = h0[1][1], lvec2.z = h0[1][2];
      lvec3.x = h0[2][0], lvec3.y = h0[2][1], lvec3.z = h0[2][2];
      set_default_recip_box();
      double vt2[3] = {vbar_eigen[0] * xdti2, vbar_eigen[1] * xdti2,
                       vbar_eigen[2] * xdti2};
      // double s1b[3] = {1.0, 1.0, 1.0};
      // double s1a[3] = {std::exp(vt2[0]), std::exp(vt2[1]), std::exp(vt2[2])};
      double s1[3] = {std::exp(vt2[0]) * sinh_id(vt2[0]),
                      std::exp(vt2[1]) * sinh_id(vt2[1]),
                      std::exp(vt2[2]) * sinh_id(vt2[2])};
      s[0] = s1[0], s[1] = s1[1], s[2] = s1[2];
   }
   for (int ir = 0; ir < nrespa; ++ir) {
      if (mid and (kw_p == KW_ATOM or kw_p == KW_MOL)) {
         double diar[3] = {std::exp(vbar_eigen[0] * xdti),
                           std::exp(vbar_eigen[1] * xdti),
                           std::exp(vbar_eigen[2] * xdti)},
                diav[3] = {s[0] * dti, s[1] * dti, s[2] * dti};
         double sr[3][3], sv[3][3];
         SymmMatrix::ODOt(sr, vbar_ortho, diar);
         SymmMatrix::ODOt(sv, vbar_ortho, diav);
         if (kw_p == KW_ATOM) {
            propagate_pos_axbv_aniso(sr, sv);
         } else if (kw_p == KW_MOL) {
            sr[0][0] = sr[0][0] - 1.0;
            sr[1][1] = sr[1][1] - 1.0;
            sr[2][2] = sr[2][2] - 1.0;
            lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
            propagate_pos_raxbv_aniso(sr, sv);
         }
      } else {
         propagate_pos(dti);
      }

      if (ir < nrespa - 1) {
         copy_pos_to_xyz(false);
         energy(vers1, RESPA_FAST, respa_tsconfig());
         if (vers1 & calc::virial) {
            lp_virial(molP);
            for (int iv = 0; iv < 9; ++iv)
               lp_vir[iv] += vir_fast[iv] / nrespa;
         }
         propagate_velocity(dti, gx, gy, gz);
      } else {
         if (constrain) {
            // TODO lp_rats1x,lp_rats1y,lp_rats1z
            lp_matvec(n, 'T', vbar_ortho, rattle_xold, rattle_yold,
                      rattle_zold);
            lp_matvec(n, 'T', vbar_ortho, xpos, ypos, zpos);
            lp_matvec(n, 'T', vbar_ortho, vx, vy, vz);
            lprat(dt, rattle_xold, rattle_yold, rattle_zold);
            lp_matvec(n, 'N', vbar_ortho, xpos, ypos, zpos);
            lp_matvec(n, 'N', vbar_ortho, vx, vy, vz);
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
   if (mid) {
      rnd_matrix[0][0] = normal<double>();
      rnd_matrix[0][1] = normal<double>();
      rnd_matrix[0][2] = normal<double>();
      rnd_matrix[1][1] = normal<double>();
      rnd_matrix[1][2] = normal<double>();
      rnd_matrix[2][2] = normal<double>();
      rnd_matrix[1][0] = rnd_matrix[0][1];
      rnd_matrix[2][0] = rnd_matrix[0][2];
      rnd_matrix[2][1] = rnd_matrix[1][2];
   }
   iso_tp_aniso(dt);
   if (constrain)
      rattle2(dt, false);


   if (molT and (istep % inform::iwrite) == 0) {
      lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
      lp_mol_kinetic();
      printf("\n"
             " Current MolKinetic     %12.4lf Kcal/mole of Frame %8d\n",
             lp_eksum, istep / inform::iwrite);
   }
}
}
