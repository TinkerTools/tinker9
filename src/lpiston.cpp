#include "box.h"
#include "energy.h"
#include "math/ou.h"
#include "md.h"
#include "random.h"
#include "rattle.h"
#include "tinker9.h"
#include "tool/error.h"
#include "tool/trimatexp.h"
#include <cmath>
#include <tinker/detail/bath.hh>
#include <tinker/detail/freeze.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void hcKinetic()
{
   auto& m = rattle_dmol;
   mdKineticEnergy(hc_eksum, hc_ekin, m.nmol, m.molmass, ratcom_vx, ratcom_vy, ratcom_vz);
}

extern void lp_mol_virial_acc();
extern void lp_mol_virial_cu();
void hcVirial()
{
   for (int iv = 0; iv < 9; ++iv)
      hc_vir[iv] = 0;
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      lp_mol_virial_cu();
   else
#endif
      lp_mol_virial_acc();
}

//====================================================================//

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

extern void hcVelIso_acc(vel_prec);
void hcVelIso(vel_prec scal)
{
   hcVelIso_acc(scal);
}

extern void hcVelAn_acc(vel_prec scal[3][3]);
void hcVelAn(vel_prec scal[3][3])
{
   hcVelAn_acc(scal);
}

void hcCenterOfMass_acc(
   const pos_prec*, const pos_prec*, const pos_prec*, pos_prec*, pos_prec*, pos_prec*);
void hcCenterOfMass(const pos_prec* atomx, const pos_prec* atomy, const pos_prec* atomz,
   pos_prec* molx, pos_prec* moly, pos_prec* molz)
{
   static_assert(
      std::is_same<pos_prec, vel_prec>::value, "pos_prec and vel_prec must be the same type.");
   hcCenterOfMass_acc(atomx, atomy, atomz, molx, moly, molz);
}

//====================================================================//

// static int nrespa, nbaro;
// static int nprtpres, iprtpres;
// static double rnd;
// static double rnd_matrix[3][3];
// static double g0;
// static double g1;
// static double D;
// static bool atomT;
// static bool molT, atomP, molP, constrain, aniso;
// static enum { KW_NULL = 0, KW_ATOM, KW_MOL } kw_p;
// static bool __pedantic;
// static const char* fmt_current_pres = "\n"
//                                       " Current Pressure       %12.4lf Atm at Step %8d\n";

void vv_lpiston_init()
{
   // auto o = stdout;

   // // RESPA keyword:    "RESPA-INNER  0.25"
   // // Barostat keyword: "BARO-OUTER    6.0"
   // nrespa = mdstuf::nrespa;
   // std::string itg;
   // get_kv("INTEGRATOR", itg, "");
   // if (itg != "RESPA")
   //    nrespa = 1;
   // get_kv("VOLUME-TRIAL", nbaro, 1);
   // nbaro = std::max(1, nbaro);

   // get_kv("FRICTION", stodyn::friction, 20.0);

   // // keyword: "PRINTOUT 5"
   // get_kv("PRINTOUT", nprtpres, 0);
   // iprtpres = 0;

   // rnd = 0.0;
   // for (int i = 0; i < 3; ++i)
   //    for (int j = 0; j < 3; ++j)
   //       rnd_matrix[i][j] = 0.0;

   // constrain = useRattle();
   // D = 3.0;

   // // keyword: "--PEDANTIC"
   // get_kbool("--PEDANTIC", __pedantic, false);
   // __pedantic = __pedantic and (not constrain);

   // std::string volscale;
   // atomP = false, molP = false;
   // atomT = false, molT = false;
   // if (constrain) {
   //    kw_p = KW_MOL;
   //    volscale = "MOLECULAR";
   //    int val = std::max(rattle_dmol.nmol, 2);
   //    g1 = D * (val - 1);
   //    g0 = D * (rattle_dmol.nmol - 1);
   //    molP = true;
   //    molT = true;
   // } else {
   //    kw_p = KW_ATOM;
   //    volscale = "ATOMIC";
   //    int val = std::max(n, 2);
   //    g1 = D * (val - 1);
   //    g0 = D * (n - 1);
   //    atomP = true;
   //    atomT = true;
   // }

   // // Nose-Hoover Chain
   // double ekt = units::gasconst * bath::kelvin;
   // vbar = 0;
   // gbar = 0;
   // qbar = (g1 + 1) * ekt * bath::taupres * bath::taupres;
   // for (int i = 0; i < maxnose; ++i) {
   //    vnh[i] = 0;
   //    gnh[i] = 0;
   //    qnh[i] = ekt * bath::tautemp * bath::tautemp;
   // }
   // qnh[0] = g0 * ekt * bath::tautemp * bath::tautemp;
   // if (nrespa > 1) {
   //    darray::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);

   //    // save fast gradients to gx1 etc.
   //    energy(calc::grad, RESPA_FAST, mdRespaTsconfig());
   //    darray::copy(g::q0, n, gx1, gx);
   //    darray::copy(g::q0, n, gy1, gy);
   //    darray::copy(g::q0, n, gz1, gz);

   //    // save slow gradients to gx2 etc.
   //    energy(calc::grad, RESPA_SLOW, mdRespaTsconfig());
   //    darray::copy(g::q0, n, gx2, gx);
   //    darray::copy(g::q0, n, gy2, gy);
   //    darray::copy(g::q0, n, gz2, gz);
   // }
   // // calculate total gradient and atomic virial
   // energy(calc::grad + calc::virial);
   // lp_virial(molP);

   // if (constrain) {
   //    darray::allocate(buffer_size(), &hc_vir_buf);
   // }

   // // isotropic NPT
   // aniso = bath::anisotrop;
   // if (aniso) {
   //    qbar /= D;
   //    vbar_matrix[0][0] = 0;
   //    vbar_matrix[0][1] = 0;
   //    vbar_matrix[0][2] = 0;
   //    vbar_matrix[1][0] = 0;
   //    vbar_matrix[1][1] = 0;
   //    vbar_matrix[1][2] = 0;
   //    vbar_matrix[2][0] = 0;
   //    vbar_matrix[2][1] = 0;
   //    vbar_matrix[2][2] = 0;
   // }

   // print(o, "\n");
   // print(o, " Friction                     %12.4lf /ps\n", stodyn::friction);
   // print(o, " Time constant for const-T    %12.4lf ps\n", bath::tautemp);
   // print(o, " Time constant for const-P    %12.4lf ps\n", bath::taupres);
   // print(o, " Pressure estimator           %12s\n", volscale);
   // print(o, " Kinetic DOF                  %12.0lf\n", g0);
   // print(o, " Pressure DOF                 %12.0lf\n", g1);
   // print(o, " NRESPA                       %12d\n", nrespa);
   // print(o, " NBARO                        %12d\n", nbaro);
   // if (aniso)
   //    print(o, " ANISOTROPIC\n");
   // if (__pedantic)
   //    print(o, " PEDANTIC\n");
   // print(o, "\n");
}

void vv_lpiston_destory()
{
   // if (constrain)
   //    darray::deallocate(hc_vir_buf);
   // if (nrespa > 1)
   //    darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);
}

//====================================================================//

// static void iso_tp(time_prec dt, bool prtpres = false, int iCurrentStep = 0)
// {
// constexpr int ns = 2;
// const double kbt = units::gasconst * bath::kelvin;
// time_prec t, t2, t4, t8, xt4;
// t = dt / ns, t2 = t / 2, t4 = t / 4, t8 = t / 8, xt4 = nbaro * t4;

// double opgxt4, omgxt4, sd;
// opgxt4 = 1.0 + stodyn::friction * xt4;
// omgxt4 = 1.0 - stodyn::friction * xt4;
// sd = std::sqrt(nbaro * dt * 2.0 * stodyn::friction * kbt / qbar) / (4 * ns);

// const virial_prec tr_vir = hc_vir[0] + hc_vir[4] + hc_vir[8];
// const double vol0 = boxVolume();
// if (atomT) {
//    lp_atom_kinetic();
// } else if (molT) {
//    lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
//    lp_mol_kinetic();
// }

// double eksum0 = hc_eksum, eksum1;
// double velsc0 = 1.0, velsc1;
// for (int k = 0; k < ns; ++k) {
//    double DelP;

//    eksum1 = eksum0;
//    velsc1 = velsc0;

//    for (int i = maxnose - 1; i > -1; --i) {
//       if (i == 0)
//          gnh[i] = (2 * eksum1 - g0 * kbt) / qnh[i];
//       else
//          gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];

//       if (i == maxnose - 1)
//          vnh[i] += gnh[i] * t4;
//       else {
//          double exptm = std::exp(-vnh[i + 1] * t8);
//          vnh[i] = (vnh[i] * exptm + gnh[i] * t4) * exptm;
//       }
//    }

//    if (atomP or molP) {
//       DelP = (1.0 + D / g1) * 2 * eksum1 - tr_vir;
//       DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
//       gbar = DelP / qbar;
//       vbar = gbar * xt4 + omgxt4 * vbar + sd * rnd;
//    }

//    double scal;
//    if (atomP or molP)
//       scal = std::exp(-t2 * (vnh[0] + (1.0 + D / g1) * vbar * nbaro));
//    else
//       scal = std::exp(-t2 * vnh[0]);
//    velsc1 *= scal;
//    eksum1 *= (scal * scal);

//    if (atomP or molP) {
//       DelP = (1.0 + D / g1) * 2 * eksum1 - tr_vir;
//       DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
//       gbar = DelP / qbar;
//       vbar = (gbar * xt4 + vbar + sd * rnd) / opgxt4;
//    }

//    for (int i = 0; i < maxnose; ++i) {
//       if (i == 0)
//          gnh[i] = (2 * eksum1 - g0 * kbt) / qnh[i];
//       else
//          gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];

//       if (i == maxnose - 1)
//          vnh[i] += gnh[i] * t4;
//       else {
//          double exptm = std::exp(-vnh[i + 1] * t8);
//          vnh[i] = (vnh[i] * exptm + gnh[i] * t4) * exptm;
//       }
//    }

//    eksum0 = eksum1;
//    velsc0 = velsc1;
// }

// const double velsc2 = velsc0 * velsc0;
// hc_eksum *= velsc2;
// for (int ii = 0; ii < 3; ++ii)
//    for (int jj = 0; jj < 3; ++jj)
//       hc_ekin[ii][jj] *= velsc2;
// if (atomT) {
//    darray::scale(g::q0, n, velsc0, vx);
//    darray::scale(g::q0, n, velsc0, vy);
//    darray::scale(g::q0, n, velsc0, vz);
// } else if (molT) {
//    hcVelIso(velsc0 - 1.0);
// }

// if (prtpres) {
//    double pres = units::prescon * (2 * hc_eksum - tr_vir) / (D * vol0);
//    print(stdout, fmt_current_pres, pres, iCurrentStep);
// }
// }

// static void lpiston_npt_iso_pedantic(int, time_prec);
// static void lpiston_npt_aniso(int, time_prec);
void vv_lpiston_npt(int istep, time_prec dt)
{
   // if (aniso) {
   //    lpiston_npt_aniso(istep, dt);
   //    return;
   // } else if (__pedantic) {
   //    lpiston_npt_iso_pedantic(istep, dt);
   //    return;
   // }

   // bool mid = (nbaro == 1) or ((istep % nbaro) == (nbaro + 1) / 2);
   // atomP = false, molP = false;
   // if (mid) {
   //    if (kw_p == KW_ATOM)
   //       atomP = true;
   //    else if (kw_p == KW_MOL)
   //       molP = true;
   // }

   // int vers1 = rc_flag & calc::vmask;
   // if ((istep % inform::iwrite) != 0)
   //    vers1 &= ~calc::energy;
   // if (not mid)
   //    vers1 &= ~calc::virial;

   // time_prec xdt = nbaro * dt, dt2 = 0.5 * dt, xdt2 = 0.5 * xdt;
   // time_prec dti = dt / nrespa, dti2 = dt2 / nrespa;
   // time_prec xdti = xdt / nrespa, xdti2 = xdt2 / nrespa;

   // iso_tp(dt);
   // if (nrespa > 1) {
   //    // gx1: fast; gx2: slow
   //    mdVel2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
   // } else {
   //    mdVel(dt2, gx, gy, gz);
   // }

   // if (constrain) {
   //    darray::copy(g::q0, n, rattle_xold, xpos);
   //    darray::copy(g::q0, n, rattle_yold, ypos);
   //    darray::copy(g::q0, n, rattle_zold, zpos);
   // }

   // virial_prec vir_fast[9] = {0};
   // energy_prec esum_f;

   // double s = 1.0;
   // if (mid) {
   //    lvec1 *= std::exp(vbar * xdt);
   //    lvec2 *= std::exp(vbar * xdt);
   //    lvec3 *= std::exp(vbar * xdt);
   //    boxSetDefaultRecip();
   //    double vt2 = vbar * xdti2;
   //    // double s1b = 1.0;
   //    // double s1a = std::exp(vt2);
   //    double s1 = std::exp(vt2) * sinhc(vt2);
   //    s = s1;
   // }
   // for (int ir = 0; ir < nrespa; ++ir) {
   //    if (mid and kw_p == KW_ATOM) {
   //       mdPosAxbv(std::exp(vbar * xdti), s * dti);
   //    } else if (mid and kw_p == KW_MOL) {
   //       lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
   //       propagate_pos_raxbv(std::exp(vbar * xdti) - 1.0, s * dti);
   //    } else {
   //       mdPos(dti);
   //    }

   //    if (ir < nrespa - 1) {
   //       mdCopyPosToXyz(false);
   //       energy(vers1, RESPA_FAST, mdRespaTsconfig());
   //       if (vers1 & calc::virial) {
   //          lp_virial(molP);
   //          for (int i = 0; i < 9; ++i)
   //             vir_fast[i] += hc_vir[i];
   //       }
   //       mdVel(dti, gx, gy, gz);
   //    } else {
   //       if (constrain) {
   //          // lp_rats1 = 1.0 / s;
   //          // lprat(dt, rattle_xold, rattle_yold, rattle_zold);
   //       }
   //       mdCopyPosToXyz(true);
   //    }
   // }

   // if (nrespa > 1) {
   //    // fast force
   //    energy(vers1, RESPA_FAST, mdRespaTsconfig());
   //    darray::copy(g::q0, n, gx1, gx);
   //    darray::copy(g::q0, n, gy1, gy);
   //    darray::copy(g::q0, n, gz1, gz);
   //    copy_energy(vers1, &esum_f);
   //    if (vers1 & calc::virial) {
   //       lp_virial(molP);
   //       for (int i = 0; i < 9; ++i)
   //          vir_fast[i] += hc_vir[i];
   //    }

   //    // slow force
   //    energy(vers1, RESPA_SLOW, mdRespaTsconfig());
   //    darray::copy(g::q0, n, gx2, gx);
   //    darray::copy(g::q0, n, gy2, gy);
   //    darray::copy(g::q0, n, gz2, gz);
   //    if (vers1 & calc::energy)
   //       esum += esum_f;
   //    if (vers1 & calc::virial) {
   //       lp_virial(molP);
   //       for (int iv = 0; iv < 9; ++iv)
   //          hc_vir[iv] += vir_fast[iv] / nrespa;
   //    }

   //    // gx1: fast; gx2: slow
   //    mdVel2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
   // } else {
   //    energy(vers1);
   //    if (vers1 & calc::virial)
   //       lp_virial(molP);
   //    mdVel(dt2, gx, gy, gz);
   // }
   // if (mid)
   //    rnd = normal<double>();
   // bool prtpres = false;
   // int iCurrentStep = 0;
   // if ((nprtpres > 0) and mid) {
   //    iprtpres++;
   //    if (iprtpres % nprtpres == 0) {
   //       prtpres = true;
   //       iCurrentStep = istep;
   //    }
   // }
   // iso_tp(dt, prtpres, iCurrentStep);
   // if (constrain)
   //    rattle2(dt, false);

   // if (molT and (istep % inform::iwrite) == 0) {
   //    lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
   //    lp_mol_kinetic();
   //    printf("\n"
   //           " Current MolKinetic     %12.4lf Kcal/mole at Frame %8d\n",
   //       hc_eksum, istep / inform::iwrite);
   // }
}

//====================================================================//

// static void iso_tp_aniso(time_prec dt, bool prtpres = false, int iCurrentStep = 0)
// {
// constexpr int ns = 2;
// const double kbt = units::gasconst * bath::kelvin;
// time_prec t, t2, t4, t8, xt4;
// t = dt / ns, t2 = t / 2, t4 = t / 4, t8 = t / 8, xt4 = nbaro * t4;

// double opgxt4, omgxt4, sd; // one plus, one minus, std dev
// opgxt4 = 1.0 + stodyn::friction * xt4;
// omgxt4 = 1.0 - stodyn::friction * xt4;
// sd = std::sqrt(nbaro * dt * 2.0 * stodyn::friction * kbt / qbar) / (4 * ns);

// const double vol0 = boxVolume();
// if (atomT) {
//    lp_atom_kinetic();
// } else if (molT) {
//    lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
//    lp_mol_kinetic();
// }

// double eksum1 = hc_eksum,
//        ekin0[3][3] = {{hc_ekin[0][0], hc_ekin[0][1], hc_ekin[0][2]},
//           {hc_ekin[1][0], hc_ekin[1][1], hc_ekin[1][2]},
//           {hc_ekin[2][0], hc_ekin[2][1], hc_ekin[2][2]}};
// double velsc0[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}, velsc1[3][3];
// for (int k = 0; k < ns; ++k) {
//    double DelP[3][3];

//    for (int i = 0; i < 3; ++i)
//       for (int j = 0; j < 3; ++j)
//          velsc1[i][j] = velsc0[i][j];

//    for (int i = maxnose - 1; i > -1; --i) {
//       if (i == 0)
//          gnh[i] = (2 * eksum1 - g0 * kbt) / qnh[i];
//       else
//          gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];

//       if (i == maxnose - 1)
//          vnh[i] += gnh[i] * t4;
//       else {
//          double exptm = std::exp(-vnh[i + 1] * t8);
//          vnh[i] = (vnh[i] * exptm + gnh[i] * t4) * exptm;
//       }
//    }

//    if (atomP or molP) {
//       for (int i = 0; i < 3; ++i) {
//          for (int j = i; j < 3; ++j) {
//             DelP[i][j] = 2 * ekin0[i][j] - hc_vir[3 * i + j];
//             if (j == i) {
//                DelP[i][i] += 2 * eksum1 / g1 - vol0 * bath::atmsph / units::prescon;
//             }
//             DelP[i][j] /= qbar;
//             vbar_matrix[i][j] =
//                DelP[i][j] * xt4 + omgxt4 * vbar_matrix[i][j] + sd * rnd_matrix[i][j];
//          }
//       }
//       if (box_shape == BoxShape::MONO) {
//          vbar_matrix[0][1] = 0, vbar_matrix[1][2] = 0;
//       } else if (box_shape == BoxShape::ORTHO or box_shape == BoxShape::OCT) {
//          vbar_matrix[0][1] = 0, vbar_matrix[1][2] = 0;
//          vbar_matrix[0][2] = 0;
//       }
//    }

//    if (atomP or molP) {
//       double tr;
//       tr = (vbar_matrix[0][0] + vbar_matrix[1][1] + vbar_matrix[2][2]) / g1;
//       // mat = -(X*vg^T + X*Tr[vg]/g1 + v1)
//       double mat[3][3] = {0};
//       for (int i = 0; i < 3; ++i) {
//          for (int j = i; j < 3; ++j) {
//             mat[i][j] = -vbar_matrix[i][j] * nbaro; // vg in use, not vg^T
//          }
//          mat[i][i] -= (tr * nbaro + vnh[0]);
//       }
//       // scal = exp(mat t/2)
//       double scal[3][3];
//       trimatExp(scal, mat, t2);
//       std::swap(scal[0][1], scal[1][0]); // transpose
//       std::swap(scal[0][2], scal[2][0]);
//       std::swap(scal[1][2], scal[2][1]);
//       matmul3(velsc1, scal);

//       // kinetic tensor = (s u).(s u)^T
//       double a = scal[0][0];
//       double b = scal[1][0], c = scal[1][1];
//       double d = scal[2][0], e = scal[2][1], f = scal[2][2];
//       double u = ekin0[0][0];
//       double v = ekin0[1][0], w = ekin0[1][1];
//       double x = ekin0[2][0], y = ekin0[2][1], z = ekin0[2][2];

//       // clang-format off
//       ekin0[0][0] = a*a*u;
//       ekin0[1][1] = b*b*u + c*c*w         + 2*b*c*v;
//       ekin0[2][2] = d*d*u + e*e*w + f*f*z + 2*d*e*v + 2*d*f*x + 2*e*f*y;
//       ekin0[1][0] = a*b*u + a*c*v;
//       ekin0[2][0] = a*d*u + a*e*v + a*f*x;
//       ekin0[2][1] = b*d*u + c*d*v + b*e*v + c*e*w + b*f*x + c*f*y;
//       ekin0[0][1] = ekin0[1][0];
//       ekin0[0][2] = ekin0[2][0];
//       ekin0[1][2] = ekin0[2][1];
//       // clang-format on
//       eksum1 = ekin0[0][0] + ekin0[1][1] + ekin0[2][2];
//    } else {
//       double scal = std::exp(-t2 * vnh[0]);
//       for (int i = 0; i < 3; ++i)
//          velsc1[i][i] *= scal;
//       scal *= scal;
//       for (int i = 0; i < 3; ++i)
//          for (int j = 0; j < 3; ++j)
//             ekin0[i][j] *= scal;
//       eksum1 *= scal;
//    }

//    if (atomP or molP) {
//       for (int i = 0; i < 3; ++i) {
//          for (int j = i; j < 3; ++j) {
//             DelP[i][j] = 2 * ekin0[i][j] - hc_vir[3 * i + j];
//             if (j == i) {
//                DelP[i][i] += 2 * eksum1 / g1 - vol0 * bath::atmsph / units::prescon;
//             }
//             DelP[i][j] /= qbar;
//             vbar_matrix[i][j] =
//                (DelP[i][j] * xt4 + vbar_matrix[i][j] + sd * rnd_matrix[i][j]) / opgxt4;
//          }
//       }
//       if (box_shape == BoxShape::MONO) {
//          vbar_matrix[0][1] = 0, vbar_matrix[1][2] = 0;
//       } else if (box_shape == BoxShape::ORTHO or box_shape == BoxShape::OCT) {
//          vbar_matrix[0][1] = 0, vbar_matrix[1][2] = 0;
//          vbar_matrix[0][2] = 0;
//       }
//    }

//    for (int i = 0; i < maxnose; ++i) {
//       if (i == 0)
//          gnh[i] = (2 * eksum1 - g0 * kbt) / qnh[i];
//       else
//          gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];

//       if (i == maxnose - 1)
//          vnh[i] += gnh[i] * t4;
//       else {
//          double exptm = std::exp(-vnh[i + 1] * t8);
//          vnh[i] = (vnh[i] * exptm + gnh[i] * t4) * exptm;
//       }
//    }

//    for (int i = 0; i < 3; ++i)
//       for (int j = 0; j < 3; ++j)
//          velsc0[i][j] = velsc1[i][j];
// }

// for (int i = 0; i < 3; ++i)
//    for (int j = 0; j < 3; ++j)
//       hc_ekin[i][j] = ekin0[i][j];
// hc_eksum = hc_ekin[0][0] + hc_ekin[1][1] + hc_ekin[2][2];
// if (atomT) {
//    lp_matvec(n, 'N', velsc0, vx, vy, vz);
// } else if (molT) {
//    vel_prec scal[3][3];
//    for (int i = 0; i < 3; ++i) {
//       for (int j = 0; j < 3; ++j) {
//          scal[i][j] = velsc0[i][j];
//       }
//       scal[i][i] -= 1.0;
//    }
//    hcVelAn(scal);
// }

// if (prtpres) {
//    double pres[3][3];
//    for (int i = 0; i < 3; ++i) {
//       for (int j = 0; j < 3; ++j) {
//          pres[i][j] = units::prescon * (2 * hc_ekin[i][j] - hc_vir[3 * i + j]) / vol0;
//       }
//    }
//    print(stdout,
//       "\n"
//       " Current Pressure X at Step %8d %12.4lf%12.4lf%12.4lf\n",
//       iCurrentStep, pres[0][0], pres[0][1], pres[0][2]);
//    print(stdout, " Current Pressure Y at Step %8d %12.4lf%12.4lf%12.4lf\n", iCurrentStep,
//       pres[1][0], pres[1][1], pres[1][2]);
//    print(stdout, " Current Pressure Z at Step %8d %12.4lf%12.4lf%12.4lf Atm\n", iCurrentStep,
//       pres[2][0], pres[2][1], pres[2][2]);
// }
// }

// static void lpiston_npt_aniso(int istep, time_prec dt)
// {
// bool mid = (nbaro == 1) or ((istep % nbaro) == (nbaro + 1) / 2);
// atomP = false, molP = false;
// if (mid) {
//    if (kw_p == KW_ATOM)
//       atomP = true;
//    else if (kw_p == KW_MOL)
//       molP = true;
// }

// int vers1 = rc_flag & calc::vmask;
// if ((istep % inform::iwrite) != 0)
//    vers1 &= ~calc::energy;
// if (not mid)
//    vers1 &= ~calc::virial;

// time_prec dt2 = 0.5 * dt, dti = dt / nrespa, dti2 = dt2 / nrespa;

// iso_tp_aniso(dt);
// if (nrespa > 1) {
//    // gx1: fast; gx2: slow
//    mdVel2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
// } else {
//    mdVel(dt2, gx, gy, gz);
// }

// if (constrain) {
//    darray::copy(g::q0, n, rattle_xold, xpos);
//    darray::copy(g::q0, n, rattle_yold, ypos);
//    darray::copy(g::q0, n, rattle_zold, zpos);
// }

// virial_prec vir_fast[9] = {0};
// energy_prec esum_f;

// if (mid) {
//    double scal[3][3];
//    trimatExp(scal, vbar_matrix, nbaro * dt);
//    double h0[3][3] = {
//       {lvec1.x, lvec1.y, lvec1.z}, {lvec2.x, lvec2.y, lvec2.z}, {lvec3.x, lvec3.y, lvec3.z}};
//    matmul3(h0, scal);
//    lvec1.x = h0[0][0], lvec1.y = h0[0][1], lvec1.z = h0[0][2];
//    lvec2.x = h0[1][0], lvec2.y = h0[1][1], lvec2.z = h0[1][2];
//    lvec3.x = h0[2][0], lvec3.y = h0[2][1], lvec3.z = h0[2][2];
//    boxSetDefaultRecip();
// }
// for (int ir = 0; ir < nrespa; ++ir) {
//    if (mid and (kw_p == KW_ATOM or kw_p == KW_MOL)) {
//       double sr[3][3], sv[3][3] = {0.0};
//       trimatExp(sr, vbar_matrix, nbaro * dti);
//       if (kw_p == KW_ATOM) {
//          double vbm[3][3];
//          for (int i = 0; i < 3; ++i)
//             for (int j = 0; j < 3; ++j)
//                vbm[i][j] = nbaro * vbar_matrix[i][j];
//          trimatTExpm1c(sv, vbm, dti);
//          propagate_pos_axbv_aniso(sr, sv);
//       } else if (kw_p == KW_MOL) {
//          sr[0][0] = sr[0][0] - 1.0;
//          sr[1][1] = sr[1][1] - 1.0;
//          sr[2][2] = sr[2][2] - 1.0;
//          // S_b(dti)*dti = I*dt/nrespa
//          sv[0][0] = dti, sv[1][1] = dti, sv[2][2] = dti;
//          lp_center_of_mass(xpos, ypos, zpos, ratcom_x, ratcom_y, ratcom_z);
//          propagate_pos_raxbv_aniso(sr, sv);
//       }
//    } else {
//       mdPos(dti);
//    }

//    if (ir < nrespa - 1) {
//       mdCopyPosToXyz(false);
//       energy(vers1, RESPA_FAST, mdRespaTsconfig());
//       if (vers1 & calc::virial) {
//          lp_virial(molP);
//          for (int iv = 0; iv < 9; ++iv)
//             hc_vir[iv] += vir_fast[iv] / nrespa;
//       }
//       mdVel(dti, gx, gy, gz);
//    } else {
//       if (constrain) {
//          // if using S_b, "rattle" will work here.
//          rattle(dt, rattle_xold, rattle_yold, rattle_zold);
//       }
//       mdCopyPosToXyz(true);
//    }
// }

// if (nrespa > 1) {
//    // fast force
//    energy(vers1, RESPA_FAST, mdRespaTsconfig());
//    darray::copy(g::q0, n, gx1, gx);
//    darray::copy(g::q0, n, gy1, gy);
//    darray::copy(g::q0, n, gz1, gz);
//    copy_energy(vers1, &esum_f);
//    if (vers1 & calc::virial) {
//       lp_virial(molP);
//       for (int i = 0; i < 9; ++i)
//          vir_fast[i] += hc_vir[i];
//    }

//    // slow force
//    energy(vers1, RESPA_SLOW, mdRespaTsconfig());
//    darray::copy(g::q0, n, gx2, gx);
//    darray::copy(g::q0, n, gy2, gy);
//    darray::copy(g::q0, n, gz2, gz);
//    if (vers1 & calc::energy)
//       esum += esum_f;
//    if (vers1 & calc::virial) {
//       lp_virial(molP);
//       for (int iv = 0; iv < 9; ++iv)
//          hc_vir[iv] += vir_fast[iv] / nrespa;
//    }

//    // gx1: fast; gx2: slow
//    mdVel2(dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
// } else {
//    energy(vers1);
//    if (vers1 & calc::virial)
//       lp_virial(molP);
//    mdVel(dt2, gx, gy, gz);
// }
// if (mid) {
//    rnd_matrix[0][0] = normal<double>();
//    rnd_matrix[0][1] = normal<double>();
//    rnd_matrix[0][2] = normal<double>();
//    rnd_matrix[1][1] = normal<double>();
//    rnd_matrix[1][2] = normal<double>();
//    rnd_matrix[2][2] = normal<double>();
//    rnd_matrix[1][0] = rnd_matrix[0][1];
//    rnd_matrix[2][0] = rnd_matrix[0][2];
//    rnd_matrix[2][1] = rnd_matrix[1][2];
// }
// bool prtpres = false;
// int iCurrentStep = 0;
// if ((nprtpres > 0) and mid) {
//    iprtpres++;
//    if (iprtpres % nprtpres == 0) {
//       prtpres = true;
//       iCurrentStep = istep;
//    }
// }
// iso_tp_aniso(dt, prtpres, iCurrentStep);
// if (constrain)
//    rattle2(dt, false);

// if (molT and (istep % inform::iwrite) == 0) {
//    lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
//    lp_mol_kinetic();
//    printf("\n"
//           " Current MolKinetic     %12.4lf Kcal/mole at Frame %8d\n",
//       hc_eksum, istep / inform::iwrite);
// }
// }

// extern void propagate_velocity_06_acc(double vbar, time_prec dt, const grad_prec* grx,
//    const grad_prec* gry, const grad_prec* grz, time_prec dt2, const grad_prec* grx2,
//    const grad_prec* gry2, const grad_prec* grz2);
// static void propagate_velocity_06(double vbar, time_prec dt, const grad_prec* grx,
//    const grad_prec* gry, const grad_prec* grz, time_prec dt2, const grad_prec* grx2,
//    const grad_prec* gry2, const grad_prec* grz2)
// {
//    propagate_velocity_06_acc(vbar, dt, grx, gry, grz, dt2, grx2, gry2, grz2);
// }
// static double* local_kinetic_atom()
// {
//    lp_atom_kinetic();
//    return &hc_eksum;
// };
// static void local_scale_vel_atom(double velsc0)
// {
//    darray::scale(g::q0, n, velsc0, vx);
//    darray::scale(g::q0, n, velsc0, vy);
//    darray::scale(g::q0, n, velsc0, vz);
// }
// separate the update for vbar in two steps
// #define TINKER_LPISTON_SEP_VBAR 0
#define TINKER_LPISTON_SEP_VBAR 1
// #define TINKER_LPISTON_VBAR_CHAIN 0
#define TINKER_LPISTON_VBAR_CHAIN 1
// double vbar_nhc[maxnose], qbar_nhc[maxnose], vbar_ekin;
// static double* local_vbar_ekin()
// {
//    vbar_ekin = 0.5 * vbar_nhc[0] * vbar_nhc[0] * qbar_nhc[0];
//    return &vbar_ekin;
// }
// static void local_scale_vbar(double velsc0)
// {
//    vbar_nhc[0] *= velsc0;
// }
// static void lpiston_npt_iso_pedantic(int istep, time_prec dt)
// {
// bool mid = (nbaro == 1) or ((istep % nbaro) == (nbaro + 1) / 2);
// atomP = false;
// if (mid) {
//    if (kw_p == KW_ATOM)
//       atomP = true;
// }

// int vers1 = rc_flag & calc::vmask;
// if ((istep % inform::iwrite) != 0)
//    vers1 &= ~calc::energy;
// if (not mid)
//    vers1 &= ~calc::virial;

// double dt2 = 0.5 * dt, dti = dt / nrespa, dti2 = dt2 / nrespa;
// double xdt = nbaro * dt, xdt2 = 0.5 * xdt;
// double xdti = xdt / nrespa, xdti2 = 0.5 * xdti;
// double b = 2 * stodyn::friction * units::gasconst * bath::kelvin / qbar;

// if (TINKER_LPISTON_VBAR_CHAIN and istep == 1) {
//    for (int i = 0; i < maxnose; ++i) {
//       vbar_nhc[i] = 0;
//       qbar_nhc[i] = qbar;
//    }
// }

// if (TINKER_LPISTON_SEP_VBAR and atomP) {
//    if (TINKER_LPISTON_VBAR_CHAIN) {
//       nhc_isot_96(xdt, maxnose, vbar_nhc, qbar_nhc, 1.0, local_vbar_ekin, local_scale_vbar);
//       vbar = vbar_nhc[0];
//    } else {
//       vbar = ornstein_uhlenbeck_process(xdt2, vbar, stodyn::friction, 0.0, b, rnd);
//    }
// }

// nhc_isot_96(dt, maxnose, vnh, qnh, g0, local_kinetic_atom, local_scale_vel_atom);

// if (atomP) {
//    virial_prec tr_vir = hc_vir[0] + hc_vir[4] + hc_vir[8];
//    double vol0 = boxVolume();
//    double DelP = (1.0 + D / g1) * 2 * hc_eksum - tr_vir;
//    DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
//    gbar = DelP / qbar;
//    if (TINKER_LPISTON_SEP_VBAR)
//       vbar = vbar + gbar * xdt2;
//    else
//       vbar = ornstein_uhlenbeck_process(xdt2, vbar, stodyn::friction, gbar, b, rnd);
// }

// if (nrespa > 1) {
//    // gx1: fast; gx2: slow
//    propagate_velocity_06((1 + D / g1) * vbar, dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
// } else {
//    propagate_velocity_06((1 + D / g1) * vbar, dt2, gx, gy, gz, 0, nullptr, nullptr, nullptr);
// }

// virial_prec vir_fast[9] = {0};
// energy_prec esum_f;

// double s = 1.0;
// if (mid) {
//    lvec1 *= std::exp(vbar * xdt);
//    lvec2 *= std::exp(vbar * xdt);
//    lvec3 *= std::exp(vbar * xdt);
//    boxSetDefaultRecip();
//    double vt2 = vbar * xdti2;
//    double s1 = std::exp(vt2) * sinhc(vt2);
//    s = s1;
// }
// for (int ir = 0; ir < nrespa; ++ir) {
//    if (mid and kw_p == KW_ATOM) {
//       mdPosAxbv(std::exp(vbar * xdti), s * dti);
//    } else {
//       mdPos(dti);
//    }

//    if (ir < nrespa - 1) {
//       mdCopyPosToXyz(false);
//       energy(vers1, RESPA_FAST, mdRespaTsconfig());
//       if (vers1 & calc::virial) {
//          lp_virial(molP);
//          for (int i = 0; i < 9; ++i)
//             vir_fast[i] += hc_vir[i];
//       }
//       mdVel(dti, gx, gy, gz);
//    } else {
//       mdCopyPosToXyz(true);
//    }
// }

// if (nrespa > 1) {
//    // fast force
//    energy(vers1, RESPA_FAST, mdRespaTsconfig());
//    darray::copy(g::q0, n, gx1, gx);
//    darray::copy(g::q0, n, gy1, gy);
//    darray::copy(g::q0, n, gz1, gz);
//    copy_energy(vers1, &esum_f);
//    if (vers1 & calc::virial) {
//       lp_virial(molP);
//       for (int i = 0; i < 9; ++i)
//          vir_fast[i] += hc_vir[i];
//    }

//    // slow force
//    energy(vers1, RESPA_SLOW, mdRespaTsconfig());
//    darray::copy(g::q0, n, gx2, gx);
//    darray::copy(g::q0, n, gy2, gy);
//    darray::copy(g::q0, n, gz2, gz);
//    if (vers1 & calc::energy)
//       esum += esum_f;
//    if (vers1 & calc::virial) {
//       lp_virial(molP);
//       for (int iv = 0; iv < 9; ++iv)
//          hc_vir[iv] += vir_fast[iv] / nrespa;
//    }

//    // gx1: fast; gx2: slow
//    propagate_velocity_06((1 + D / g1) * vbar, dti2, gx1, gy1, gz1, dt2, gx2, gy2, gz2);
// } else {
//    energy(vers1);
//    if (vers1 & calc::virial)
//       lp_virial(molP);
//    propagate_velocity_06((1 + D / g1) * vbar, dt2, gx, gy, gz, 0, nullptr, nullptr, nullptr);
// }

// if (mid)
//    rnd = normal<double>();
// bool prtpres = false;
// if ((nprtpres > 0) and mid) {
//    iprtpres++;
//    if (iprtpres % nprtpres == 0) {
//       prtpres = true;
//    }
// }
// if (prtpres) {
//    virial_prec tr_vir = hc_vir[0] + hc_vir[4] + hc_vir[8];
//    double vol0 = boxVolume();
//    double pres = units::prescon * (2 * hc_eksum - tr_vir) / (D * vol0);
//    print(stdout, fmt_current_pres, pres, istep);
// }

// if (atomP) {
//    virial_prec tr_vir = hc_vir[0] + hc_vir[4] + hc_vir[8];
//    double vol0 = boxVolume();
//    double DelP = (1.0 + D / g1) * 2 * hc_eksum - tr_vir;
//    DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
//    gbar = DelP / qbar;
//    if (TINKER_LPISTON_SEP_VBAR)
//       vbar = vbar + gbar * xdt2;
//    else
//       vbar = ornstein_uhlenbeck_process(xdt2, vbar, stodyn::friction, gbar, b, rnd);
// }

// nhc_isot_96(dt, maxnose, vnh, qnh, g0, local_kinetic_atom, local_scale_vel_atom);

// if (TINKER_LPISTON_SEP_VBAR and atomP) {
//    if (TINKER_LPISTON_VBAR_CHAIN) {
//       nhc_isot_96(xdt, maxnose, vbar_nhc, qbar_nhc, 1.0, local_vbar_ekin, local_scale_vbar);
//       vbar = vbar_nhc[0];
//    } else {
//       vbar = ornstein_uhlenbeck_process(xdt2, vbar, stodyn::friction, 0.0, b, rnd);
//    }
// }

// if (molT and (istep % inform::iwrite) == 0) {
//    lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
//    lp_mol_kinetic();
//    printf("\n"
//           " Current MolKinetic     %12.4lf Kcal/mole at Frame %8d\n",
//       hc_eksum, istep / inform::iwrite);
// }
// }
}
