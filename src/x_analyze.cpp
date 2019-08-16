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
  if (use_potent(polar_term))
    print(out, fmt, "Polarization", get_energy(ep), get_count(nep));

  finish();
}
TINKER_NAMESPACE_END
