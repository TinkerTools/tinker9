#include "energy.h"
#include "md.h"
#include "nblist.h"
#include "osrw.h"
#include "potent.h"
#include "tinker9.h"
#include "tool/io.h"
#include <fstream>
#include <tinker/detail/bound.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/files.hh>
#include <tinker/detail/moment.hh>
#include <tinker/detail/units.hh>

namespace tinker {
static void x_analyze_e()
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
      print(out, fmt, "Bond Stretching", energy_eb, count_bonded_term(bond_term));

   if (use_potent(angle_term))
      print(out, fmt, "Angle Bending", energy_ea, count_bonded_term(angle_term));

   if (use_potent(strbnd_term))
      print(out, fmt, "Stretch-Bend", energy_eba, count_bonded_term(strbnd_term));

   if (use_potent(urey_term))
      print(out, fmt, "Urey-Bradley", energy_eub, count_bonded_term(urey_term));

   if (use_potent(opbend_term))
      print(out, fmt, "Out-of-Plane Bend", energy_eopb, count_bonded_term(opbend_term));

   if (use_potent(improp_term))
      print(out, fmt, "Improper Dihedral", energy_eid, count_bonded_term(improp_term));

   if (use_potent(imptors_term))
      print(out, fmt, "Improper Torsion", energy_eit, count_bonded_term(imptors_term));

   if (use_potent(torsion_term))
      print(out, fmt, "Torsional Angle", energy_et, count_bonded_term(torsion_term));

   if (use_potent(pitors_term))
      print(out, fmt, "Pi-Orbital Torsion", energy_ept, count_bonded_term(pitors_term));

   if (use_potent(strtor_term))
      print(out, fmt, "Stretch-Torsion", energy_ebt, count_bonded_term(strtor_term));

   if (use_potent(angtor_term))
      print(out, fmt, "Angle-Torsion", energy_eat, count_bonded_term(angle_term));

   if (use_potent(tortor_term))
      print(out, fmt, "Torsion-Torsion", energy_ett, count_bonded_term(tortor_term));

   if (use_potent(vdw_term))
      print(out, fmt, "Van der Waals", energy_ev, count_reduce(nev));

   if (use_potent(repuls_term))
      print(out, fmt, "Repulsion", energy_er, count_reduce(nrep));

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
      print(out, fmt, "Geometric Restraints", energy_eg, count_bonded_term(geom_term));
}

void moments();
static void x_analyze_m()
{
   auto out = stdout;
   moments();
   print(out,
      "\n"
      " Total Electric Charge :%12s%13.5lf Electrons\n",
      "", moment::netchg);
   print(out,
      "\n"
      " Dipole Moment Magnitude :%10s%13.3lf Debye\n"
      "\n"
      " Dipole X,Y,Z-Components :%10s%13.3lf%13.3lf%13.3lf\n",
      "", moment::netdpl, "", moment::xdpl, moment::ydpl, moment::zdpl);
   print(out,
      "\n"
      " Quadrupole Moment Tensor :%9s%13.3lf%13.3lf%13.3lf\n"
      "      (Buckinghams)%17s%13.3lf%13.3lf%13.3lf\n"
      "%36s%13.3lf%13.3lf%13.3lf\n",
      "", moment::xxqpl, moment::xyqpl, moment::xzqpl, "", moment::yxqpl, moment::yyqpl,
      moment::yzqpl, "", moment::zxqpl, moment::zyqpl, moment::zzqpl);
   print(out,
      "\n"
      " Principal Axes Quadrupole :%8s%13.3lf%13.3lf%13.3lf\n",
      "", moment::netqpl[0], moment::netqpl[1], moment::netqpl[2]);

   if (chgpot::dielec != 1) {
      print(out,
         "\n"
         " Dielectric Constant :%14s%13.3lf\n",
         "", chgpot::dielec);
      print(out, " Effective Total Charge :%11s%13.5lf Electrons\n", "",
         moment::netchg / std::sqrt(chgpot::dielec));
      print(out, " Effective Dipole Moment :%10s%13.3lf Debye\n", "",
         moment::netdpl / std::sqrt(chgpot::dielec));
   }

   // radius of gyration and moments of inertia
   double rg;
   tinker_f_gyrate(&rg);
   print(out,
      "\n"
      " Radius of Gyration :%15s%13.3lf Angstroms\n",
      "", rg);
   int one = 1;
   tinker_f_inertia(&one);
}

static void x_analyze_v()
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

   double pres = 0;
   int temp = 298;
   const char* fmt_p = "\n"
                       " Pressure (Temp %3d K) :            %13.3lf Atmospheres\n";
   if (bound::use_bounds) {
      double vol = boxVolume();
      double tr_vir = vir[0] + vir[4] + vir[8];
      double pres_vir = -tr_vir;
      pres = 3 * n * units::gasconst * temp + pres_vir;
      pres *= units::prescon / (3 * vol);
      pres_vir *= units::prescon / (3 * vol);
      print(out, fmt_p, temp, pres);
      print(out, " Pressure From Virial               %13.3lf Atmospheres\n", pres_vir);
   } else {
      print(out, fmt_p, temp, pres);
   }
}

void x_analyze(int, char**)
{
   initial();
   tinker_f_getxyz();
   tinker_f_mechanic();
   mechanic2();

   char string[240];
   int exist = false;
   std::string opt;
   nextarg(string, exist);
   if (exist) {
      ioReadString(opt, string);
   }
   std::string prompt = R"(
 The Tinker Energy Analysis Utility Can :

 Total Potential Energy and its Components [E]
 Electrostatic Moments and Principle Axes [M]
 Internal Virial, dE/dV Values & Pressure [V]

 Enter the Desired Analysis Types [E] :  )";
   ioReadStream(opt, prompt, std::string("#"), [](std::string s) {
      Text::upcase(s);
      auto failed = std::string::npos;
      if (s.find("E") != failed or s.find("M") != failed or s.find("V") != failed)
         return 0;
      else
         return 1;
   });
   Text::upcase(opt);

   int flags = calc::xyz + calc::mass;
   flags += (calc::energy + calc::grad + calc::virial + calc::analyz);
   rc_flag = flags;
   initialize();

   auto failed = std::string::npos;
   auto out = stdout;
   FstrView fsw = files::filename;
   std::string fname = fsw.trim();
   std::ifstream ipt(fname);
   int done = false;
   int nframe_processed = 0;
   do {
      mdReadFrameCopyinToXyz(ipt, done);
      refresh_neighbors();
      nframe_processed++;
      if (nframe_processed > 1)
         print(out, "\n Analysis for Archive Structure :%16d\n", nframe_processed);
      if (opt.find("E") != failed)
         x_analyze_e();
      if (opt.find("M") != failed)
         x_analyze_m();
      if (opt.find("V") != failed)
         x_analyze_v();
   } while (not done);

   finish();
   tinker_f_final();
}
}
