#include "tinker_rt.h"
#include <chrono>
#include <tinker/detail/bath.hh>
#include <tinker/detail/bound.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/keys.hh>
#include <tinker/detail/mdstuf.hh>

TINKER_NAMESPACE_BEGIN
void x_dynamic(int, char**)
{
   char string[240];
   int exist = false;

   initial();
   TINKER_RT(getxyz)();
   mechanic();

   bath::kelvin = 0;
   bath::atmsph = 0;
   bath::isothermal = false;
   bath::isobaric = false;

   // check for keywords containing any altered parameters

   fstr_view integrate = mdstuf::integrate;
   integrate = "BEEMAN";
   for (auto i = 0; i < keys::nkey; ++i) {
      auto record = Text::string(keys::keyline[i]);
      auto vs = Text::split(record);
      if (vs.size()) {
         auto keyword = vs.at(0);
         Text::upcase(keyword);
         if (keyword == "INTEGRATOR") {
            record = vs.at(1);
            Text::upcase(record);
            integrate = record;
         }
      }
   }

   // initialize the simulation length as number of time steps

   int nstep = -1;
   nextarg(string, exist);
   if (exist) {
      read_string(nstep, string);
   }
   read_stream(nstep,
               "\n Enter the Number of Dynamics Steps to be Taken :  ", 0,
               [](int i) { return i < 0; });

   // get the length of the dynamics time step in picoseconds

   double dt = -1;
   nextarg(string, exist);
   if (exist) {
      read_string(dt, string);
   }
   read_stream(dt,
               "\n Enter the Time Step Length in Femtoseconds [1.0] :  ", 1.0,
               [](double i) { return i <= 0; });
   dt *= 0.001;

   // enforce bounds on thermostat and barostat coupling times

   bath::tautemp = std::max(bath::tautemp, dt);
   bath::taupres = std::max(bath::taupres, dt);

   // set the time between trajectory snapshot coordinate saves

   double dtsave = -1;
   nextarg(string, exist);
   if (exist) {
      read_string(dtsave, string);
   }
   read_stream(dtsave,
               "\n Enter Time between Saves in Picoseconds [0.1] :  ", 0.1,
               [](double i) { return i <= 0; });
   inform::iwrite = std::round(dtsave / dt);

   // get choice of statistical ensemble for periodic system

   const char* ask_kelvin =
      "\n Enter the Desired Temperature in Degrees K [298] :  ";
   const double default_kelvin = 298.0;
   auto invalid_kelvin = [](double t) { return t < 0; };
   const char* ask_atm = "\n Enter the Desired Pressure in Atm [1.0] :  ";
   const double default_atm = 1.0;
   auto invalid_atm = [](double t) { return t <= 0; };

   int mode = -1;
   if (bound::use_bounds) {
      nextarg(string, exist);
      if (exist) {
         read_string(mode, string);
      }
      std::string prompt = R"(
 Available Statistical Mechanical Ensembles :
    (1) Microcanonical (NVE)
    (2) Canonical (NVT)
    (3) Isoenthalpic-Isobaric (NPH)
    (4) Isothermal-Isobaric (NPT)
 Enter the Number of the Desired Choice  [1] :  )";
      read_stream(mode, prompt, 1, [](int i) { return i <= 0; });

      if (integrate == "BUSSI" || integrate == "NOSE-HOOVER" ||
          integrate == "GHMC") {
         if (mode != 4) {
            mode = 4;
            print(
               stdout,
               "\n Switching to NPT Ensemble as Required by Chosen Integrator");
         }
      }

      if (mode == 2 || mode == 4) {
         bath::isothermal = true;
         bath::kelvin = -1;
         nextarg(string, exist);
         if (exist) {
            read_string(bath::kelvin, string);
         }
         read_stream(bath::kelvin, ask_kelvin, default_kelvin, invalid_kelvin);
      }

      if (mode == 3 || mode == 4) {
         bath::isobaric = true;
         bath::atmsph = -1;
         nextarg(string, exist);
         if (exist) {
            read_string(bath::atmsph, string);
         }
         read_stream(bath::atmsph, ask_atm, default_atm, invalid_atm);
      }
   } else {
      nextarg(string, exist);
      if (exist) {
         read_string(mode, string);
      }
      std::string prompt = R"(
 Available Simulation Control Modes :
    (1) Constant Total Energy Value (E)
    (2) Constant Temperature via Thermostat (T)
 Enter the Number of the Desired Choice [1] :  )";
      read_stream(mode, prompt, 1, [](int i) { return i <= 0; });

      if (mode == 2) {
         bath::isothermal = true;
         bath::kelvin = -1;
         nextarg(string, exist);
         if (exist) {
            read_string(bath::kelvin, string);
         }
         read_stream(bath::kelvin, ask_kelvin, default_kelvin, invalid_kelvin);
      }
   }

   // perform the setup functions needed to run dynamics

   TINKER_RT(mdinit)();

   int flags = calc::md;
   flags += (calc::xyz + calc::vel + calc::mass);
   flags += (calc::energy + calc::grad);
   if (bath::isobaric == true)
      flags += calc::virial;

   rc_flag = flags;
   initialize();

   if (integrate == "VERLET")
      print(stdout,
            " Molecular Dynamics Trajectory via"
            " Velocity Verlet Algorithm\n");
   else if (integrate == "RESPA")
      print(stdout,
            " Molecular Dynamics Trajectory via"
            " r-RESPA MTS Algorithm\n");

   auto t_start = std::chrono::steady_clock::now();
   propagate(nstep, dt);
   auto t_end = std::chrono::steady_clock::now();

   auto d_us =
      std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start)
         .count();
   double us1 = (dt * 1000.) * nstep * 86400;

   const char* fmt_flt = " %-14s%-9s%18.4f\n";
   const char* fmt_int = " %-14s%-9s%18d\n";
   print(stdout, "\n");
   print(stdout, fmt_flt, "Performance:", "ns/day", us1 / d_us);
   print(stdout, fmt_flt, "", "Wall Time", d_us / 1000000.);
   print(stdout, fmt_int, "", "Steps", nstep);
   print(stdout, fmt_int, "", "Updates", nstep / inform::iwrite);
   print(stdout, fmt_flt, "", "Time Step", dt * 1000);
   print(stdout, fmt_int, "", "Atoms", n);

   finish();

   // perform any final tasks before program exit

   TINKER_RT(final)();
}
TINKER_NAMESPACE_END
