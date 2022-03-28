#include "ff/amoeba/elecamoeba.h"
#include "ff/atom.h"
#include "ff/box.h"
#include "ff/energy.h"
#include "ff/hippo/edisp.h"
#include "ff/hippo/elechippo.h"
#include "ff/hippo/erepel.h"
#include "ff/nblist.h"
#include "ff/pchg/echarge.h"
#include "ff/pchg/evalence.h"
#include "ff/pchg/evdw.h"
#include "ff/potent.h"
#include "md/inc.h"
#include "md/osrw.h"
#include "tinker9.h"
#include "tool/darray.h"
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

   if (usePotent(Potent::BOND))
      print(out, fmt, "Bond Stretching", energy_eb, countBondedTerm(Potent::BOND));

   if (usePotent(Potent::ANGLE))
      print(out, fmt, "Angle Bending", energy_ea, countBondedTerm(Potent::ANGLE));

   if (usePotent(Potent::STRBND))
      print(out, fmt, "Stretch-Bend", energy_eba, countBondedTerm(Potent::STRBND));

   if (usePotent(Potent::UREY))
      print(out, fmt, "Urey-Bradley", energy_eub, countBondedTerm(Potent::UREY));

   if (usePotent(Potent::OPBEND))
      print(out, fmt, "Out-of-Plane Bend", energy_eopb, countBondedTerm(Potent::OPBEND));

   if (usePotent(Potent::IMPROP))
      print(out, fmt, "Improper Dihedral", energy_eid, countBondedTerm(Potent::IMPROP));

   if (usePotent(Potent::IMPTORS))
      print(out, fmt, "Improper Torsion", energy_eit, countBondedTerm(Potent::IMPTORS));

   if (usePotent(Potent::TORSION))
      print(out, fmt, "Torsional Angle", energy_et, countBondedTerm(Potent::TORSION));

   if (usePotent(Potent::PITORS))
      print(out, fmt, "Pi-Orbital Torsion", energy_ept, countBondedTerm(Potent::PITORS));

   if (usePotent(Potent::STRTOR))
      print(out, fmt, "Stretch-Torsion", energy_ebt, countBondedTerm(Potent::STRTOR));

   if (usePotent(Potent::ANGTOR))
      print(out, fmt, "Angle-Torsion", energy_eat, countBondedTerm(Potent::ANGLE));

   if (usePotent(Potent::TORTOR))
      print(out, fmt, "Torsion-Torsion", energy_ett, countBondedTerm(Potent::TORTOR));

   if (usePotent(Potent::VDW))
      print(out, fmt, "Van der Waals", energy_ev, countReduce(nev));

   if (usePotent(Potent::REPULS))
      print(out, fmt, "Repulsion", energy_er, countReduce(nrep));

   if (usePotent(Potent::DISP))
      print(out, fmt, "Dispersion", energy_edsp, countReduce(ndisp));

   if (usePotent(Potent::CHARGE))
      print(out, fmt, "Charge-Charge", energy_ec, countReduce(nec));

   if (usePotent(Potent::MPOLE))
      print(out, fmt, "Atomic Multipoles", energy_em, countReduce(nem));

   if (usePotent(Potent::POLAR))
      print(out, fmt, "Polarization", energy_ep, countReduce(nep));

   if (usePotent(Potent::CHGTRN))
      print(out, fmt, "Charge Transfer", energy_ect, countReduce(nct));

   if (usePotent(Potent::GEOM))
      print(out, fmt, "Geometric Restraints", energy_eg, countBondedTerm(Potent::GEOM));
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

void xAnalyze(int, char**)
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
      nblistRefresh();
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
