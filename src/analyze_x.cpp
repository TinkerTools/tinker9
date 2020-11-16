#include "energy.h"
#include "md.h"
#include "nblist.h"
#include "osrw.h"
#include "potent.h"
#include "tinker_rt.h"
#include "tool/io_print.h"
#include "tool/io_read.h"
#include <fstream>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/files.hh>
#include <tinker/detail/mpole.hh>
#include <tinker/detail/polar.hh>

namespace tinker {
void x_analyze_e();
void moments();
void x_analyze_m();
void x_analyze_v();

void x_analyze(int, char**)
{
   initial();
   TINKER_RT(getxyz)();
   mechanic();
   mechanic2();

   char string[240];
   int exist = false;
   std::string opt;
   nextarg(string, exist);
   if (exist) {
      read_string(opt, string);
   }
   //   std::string prompt = R"(
   //  The Tinker Energy Analysis Utility Can :

   //  General System and Force Field Information [G]
   //  Force Field Parameters for Interactions [P]
   //  Total Potential Energy and its Components [E]
   //  Energy Breakdown over Each of the Atoms [A]
   //  List of the Large Individual Interactions [L]
   //  Details for All Individual Interactions [D]
   //  Electrostatic Moments and Principle Axes [M]
   //  Internal Virial, dE/dV Values & Pressure [V]
   //  Connectivity Lists for Each of the Atoms [C]

   //  Enter the Desired Analysis Types [G,P,E,A,L,D,M,V,C] :  )";
   std::string prompt = R"(
 The Tinker Energy Analysis Utility Can :

 Total Potential Energy and its Components [E]
 Electrostatic Moments and Principle Axes [M]
 Internal Virial, dE/dV Values & Pressure [V]

 Enter the Desired Analysis Types [E] :  )";
   read_stream(opt, prompt, std::string("#"), [](std::string s) {
      Text::upcase(s);
      auto failed = std::string::npos;
      if (s.find("E") != failed or s.find("M") != failed or
          s.find("V") != failed)
         return 0;
      else
         return 1;
   });


   int flags = calc::xyz + calc::mass;
   flags += (calc::energy + calc::grad + calc::virial + calc::analyz);
   rc_flag = flags;
   initialize();


   auto failed = std::string::npos;
   Text::upcase(opt);


   auto out = stdout;
   fstr_view fsw = files::filename;
   std::string fname = fsw.trim();
   std::ifstream ipt(fname);
   int done = false;
   read_frame_copyin_to_xyz(ipt, done);
   refresh_neighbors();
   int nframe_processed = 1;
   if (opt.find("E") != failed)
      x_analyze_e();
   if (opt.find("M") != failed)
      x_analyze_m();
   if (opt.find("V") != failed)
      x_analyze_v();
   while (!done) {
      read_frame_copyin_to_xyz(ipt, done);
      refresh_neighbors();
      nframe_processed++;
      print(out, "\n Analysis for Archive Structure :%16d\n", nframe_processed);
      if (opt.find("E") != failed)
         x_analyze_e();
      if (opt.find("M") != failed)
         x_analyze_m();
      if (opt.find("V") != failed)
         x_analyze_v();
   }


   finish();
   TINKER_RT(final)();
}

void x_analyze_e()
{
   if (use_osrw)
      osrw_energy(calc::energy + calc::analyz);
   else
      energy(calc::energy + calc::analyz);

   auto& out = stdout;
   print(out, "\n Total Potential Energy :        %16.4f Kcal/mole\n", esum);
   print(out,
         "\n Energy Component Breakdown :           Kcal/mole        "
         "Interactions\n\n");

   const char* fmt = " %-29s %18.4f %16d\n";

   if (use_potent(bond_term))
      print(out, fmt, "Bond Stretching", energy_eb,
            count_bonded_term(bond_term));

   if (use_potent(angle_term))
      print(out, fmt, "Angle Bending", energy_ea,
            count_bonded_term(angle_term));

   if (use_potent(strbnd_term))
      print(out, fmt, "Stretch-Bend", energy_eba,
            count_bonded_term(strbnd_term));

   if (use_potent(urey_term))
      print(out, fmt, "Urey-Bradley", energy_eub, count_bonded_term(urey_term));

   if (use_potent(opbend_term))
      print(out, fmt, "Out-of-Plane Bend", energy_eopb,
            count_bonded_term(opbend_term));

   if (use_potent(imptors_term))
      print(out, fmt, "Improper Torsion", energy_eit,
            count_bonded_term(imptors_term));

   if (use_potent(torsion_term))
      print(out, fmt, "Torsional Angle", energy_et,
            count_bonded_term(torsion_term));

   if (use_potent(pitors_term))
      print(out, fmt, "Pi-Orbital Torsion", energy_ept,
            count_bonded_term(pitors_term));

   if (use_potent(tortor_term))
      print(out, fmt, "Torsion-Torsion", energy_ett,
            count_bonded_term(tortor_term));

   if (use_potent(vdw_term))
      print(out, fmt, "Van der Waals", energy_ev, count_reduce(nev));

   if (use_potent(disp_term))
      print(out, fmt, "Dispersion", energy_edsp, count_reduce(ndisp));

   if (use_potent(charge_term))
      print(out, fmt, "Charge-Charge", energy_ec, count_reduce(nec));

   if (use_potent(mpole_term))
      print(out, fmt, "Atomic Multipoles", energy_em, count_reduce(nem));

   if (use_potent(polar_term))
      print(out, fmt, "Polarization", energy_ep, count_reduce(nep));

   if (use_potent(chgtrn_term))
      print(out, fmt, "Charge Transfer", energy_ect, count_reduce(nct));

   if (use_potent(geom_term))
      print(out, fmt, "Geometric Restraints", energy_eg,
            count_bonded_term(geom_term));
}


void x_analyze_m()
{
   bounds();
   // download x y z
   std::vector<pos_prec> xv(n), yv(n), zv(n);
   darray::copyout(g::q0, n, xv.data(), xpos);
   darray::copyout(g::q0, n, yv.data(), ypos);
   darray::copyout(g::q0, n, zv.data(), zpos);
   wait_for(g::q0);
   for (int i = 0; i < n; ++i) {
      atoms::x[i] = xv[i];
      atoms::y[i] = yv[i];
      atoms::z[i] = zv[i];
   }

   // download rpole, uind
   if (use_potent(mpole_term) or use_potent(polar_term)) {
      std::vector<real> rpolev(n * 10), uindv(n * 3);
      mpole_init(calc::energy);
      darray::copyout(g::q0, n * 10, rpolev.data(), &rpole[0][0]);
      if (use_potent(polar_term)) {
         induce(uind, uinp);
         darray::copyout(g::q0, n * 3, uindv.data(), &uind[0][0]);
      }
      wait_for(g::q0);
      for (int i = 0; i < n; ++i) {
         int t = 0;
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_0];
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_x];
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_y];
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_z];
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_xx]; // xx
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_xy]; // xy
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_xz]; // xz
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_yx]; // yx
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_yy]; // yy
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_yz]; // yz
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_zx]; // zx
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_zy]; // zy
         mpole::rpole[13 * i + (t++)] = rpolev[10 * i + mpl_pme_zz]; // zz
         polar::uind[3 * i + 0] = uindv[3 * i + 0];
         polar::uind[3 * i + 1] = uindv[3 * i + 1];
         polar::uind[3 * i + 2] = uindv[3 * i + 2];
      }
   }
   moments();
}


void x_analyze_v()
{
   if (use_osrw)
      osrw_energy(calc::grad + calc::virial);
   else
      energy(calc::grad + calc::virial);
   auto& out = stdout;

   const char* fmt = " %-36s%12.3f %12.3f %12.3f\n";
   print(out, "\n");
   print(out, fmt, "Internal Virial Tensor :", vir[0], vir[1], vir[2]);
   print(out, fmt, "", vir[3], vir[4], vir[5]);
   print(out, fmt, "", vir[6], vir[7], vir[8]);
}
}
