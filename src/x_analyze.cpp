#include "energy.h"
#include "io_print.h"
#include "io_read.h"
#include "md.h"
#include "potent.h"
#include "tinker_rt.h"

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
   if (opt.find("E") != failed)
      x_analyze_e();
   if (opt.find("V") != failed)
      x_analyze_v();


   finish();
   TINKER_RT(final)();
}

void x_analyze_e()
{
   energy_potential(calc::energy + calc::analyz);

   auto& out = stdout;
   print(out, "\n Total Potential Energy :        {:16.4f} Kcal/mole\n", esum);
   print(out,
         "\n Energy Component Breakdown :           Kcal/mole        "
         "Interactions\n\n");

   const char* fmt = " {:<29s}{:19.4f}{:17d}\n";

   if (use_potent(bond_term))
      print(out, fmt, "Bond Stretching", get_energy(eb),
            count_bonded_term(bond_term));

   if (use_potent(angle_term))
      print(out, fmt, "Angle Bending", get_energy(ea),
            count_bonded_term(angle_term));

   if (use_potent(strbnd_term))
      print(out, fmt, "Stretch-Bend", get_energy(eba),
            count_bonded_term(strbnd_term));

   if (use_potent(urey_term))
      print(out, fmt, "Urey-Bradley", get_energy(eub),
            count_bonded_term(urey_term));

   if (use_potent(opbend_term))
      print(out, fmt, "Out-of-Plane Bend", get_energy(eopb),
            count_bonded_term(opbend_term));

   if (use_potent(torsion_term))
      print(out, fmt, "Torsional Angle", get_energy(et),
            count_bonded_term(torsion_term));

   if (use_potent(pitors_term))
      print(out, fmt, "Pi-Orbital Torsion", get_energy(ept),
            count_bonded_term(pitors_term));

   if (use_potent(tortor_term))
      print(out, fmt, "Torsion-Torsion", get_energy(ett),
            count_bonded_term(tortor_term));

   if (use_potent(vdw_term))
      print(out, fmt, "Van der Waals", get_energy(ev), get_count(nev));

   if (use_potent(mpole_term))
      print(out, fmt, "Atomic Multipoles", get_energy(em), get_count(nem));

   if (use_potent(polar_term))
      print(out, fmt, "Polarization", get_energy(ep), get_count(nep));
}


void x_analyze_v()
{
   energy_potential(calc::grad + calc::virial);
   auto& out = stdout;

   const char* fmt = " {:36s} {:12.3f} {:12.3f} {:12.3f}\n";
   print(out, "\n");
   print(out, fmt, "Internal Virial Tensor :", vir[0], vir[1], vir[2]);
   print(out, fmt, "", vir[3], vir[4], vir[5]);
   print(out, fmt, "", vir[6], vir[7], vir[8]);
}
TINKER_NAMESPACE_END
