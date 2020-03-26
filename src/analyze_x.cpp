#include "energy.h"
#include "io_print.h"
#include "io_read.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tinker_rt.h"
#include <fstream>
#include <tinker/detail/files.hh>

TINKER_NAMESPACE_BEGIN
void x_analyze_e();
void x_analyze_v();

void x_analyze(int argc, char** argv)
{
   TINKER_RT(initial)();
   TINKER_RT(getxyz)();
   TINKER_RT(mechanic)();

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
 Internal Virial, dE/dV Values & Pressure [V]

 Enter the Desired Analysis Types [E] :  )";
   read_stream(opt, prompt, std::string("#"), [](std::string s) {
      Text::upcase(s);
      auto failed = std::string::npos;
      if (s.find("E") != failed || s.find("V") != failed)
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
   if (opt.find("V") != failed)
      x_analyze_v();
   while (!done) {
      read_frame_copyin_to_xyz(ipt, done);
      refresh_neighbors();
      nframe_processed++;
      print(out, "\n Analysis for Archive Structure :%16d\n", nframe_processed);
      if (opt.find("E") != failed)
         x_analyze_e();
      if (opt.find("V") != failed)
         x_analyze_v();
   }


   finish();
   TINKER_RT(final)();
}

void x_analyze_e()
{
   energy(calc::energy + calc::analyz);

   auto& out = stdout;
   print(out, "\n Total Potential Energy :        %16.4 Kcal/mole\n", esum);
   print(out,
         "\n Energy Component Breakdown :           Kcal/mole        "
         "Interactions\n\n");

   const char* fmt = " %-29s %18.4f %16d\n";

   if (use_potent(bond_term))
      print(out, fmt, "Bond Stretching", energy_reduce(eb),
            count_bonded_term(bond_term));

   if (use_potent(angle_term))
      print(out, fmt, "Angle Bending", energy_reduce(ea),
            count_bonded_term(angle_term));

   if (use_potent(strbnd_term))
      print(out, fmt, "Stretch-Bend", energy_reduce(eba),
            count_bonded_term(strbnd_term));

   if (use_potent(urey_term))
      print(out, fmt, "Urey-Bradley", energy_reduce(eub),
            count_bonded_term(urey_term));

   if (use_potent(opbend_term))
      print(out, fmt, "Out-of-Plane Bend", energy_reduce(eopb),
            count_bonded_term(opbend_term));

   if (use_potent(torsion_term))
      print(out, fmt, "Torsional Angle", energy_reduce(et),
            count_bonded_term(torsion_term));

   if (use_potent(pitors_term))
      print(out, fmt, "Pi-Orbital Torsion", energy_reduce(ept),
            count_bonded_term(pitors_term));

   if (use_potent(tortor_term))
      print(out, fmt, "Torsion-Torsion", energy_reduce(ett),
            count_bonded_term(tortor_term));

   if (use_potent(vdw_term))
      print(out, fmt, "Van der Waals", energy_reduce(ev), count_reduce(nev));

   if (use_potent(charge_term))
      print(out, fmt, "Charge-Charge", energy_reduce(ec), count_reduce(nec));

   if (use_potent(mpole_term))
      print(out, fmt, "Atomic Multipoles", energy_reduce(em),
            count_reduce(nem));

   if (use_potent(polar_term))
      print(out, fmt, "Polarization", energy_reduce(ep), count_reduce(nep));

   if (use_potent(geom_term))
      print(out, fmt, "Geometric Restraints", energy_reduce(eg),
            count_bonded_term(geom_term));
}


void x_analyze_v()
{
   energy(calc::grad + calc::virial);
   auto& out = stdout;

   const char* fmt = " %-36s%12.3f %12.3f %12.3f\n";
   print(out, "\n");
   print(out, fmt, "Internal Virial Tensor :", vir[0], vir[1], vir[2]);
   print(out, fmt, "", vir[3], vir[4], vir[5]);
   print(out, fmt, "", vir[6], vir[7], vir[8]);
}
TINKER_NAMESPACE_END
