#include "energy.h"
#include "io_print.h"
#include "io_read.h"
#include "md.h"
#include "potent.h"
#include "rc_man.h"
#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
static void option_e();

void x_analyze(int argc, char** argv) {
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

 Enter the Desired Analysis Types [E] :  )";
  read_stream(opt, prompt, std::string("#"), [](std::string s) {
    Text::upcase(s);
    if (s == "E")
      return false;
    else
      return true;
  });

  Text::upcase(opt);
  if (opt == "E")
    option_e();

  TINKER_RT(final)();
}

static void option_e() {
  int flags = calc::xyz + calc::mass;
  flags += (calc::energy + calc::analyz);

  rc_flag = flags;
  initialize();

  energy_potential(rc_flag & calc::vmask);

  double sum = get_energy(esum);
  auto& out = std::cout;
  print(out, "\n Total Potential Energy :        {:16.4f} Kcal/mole\n", sum);
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

  finish();
}
TINKER_NAMESPACE_END
