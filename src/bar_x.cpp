#include "energy.h"
#include "tinker_rt.h"
#include <array>
#include <tinker/detail/files.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/keys.hh>
#include <tinker/detail/titles.hh>


namespace tinker {
namespace {
constexpr int maxframe = 1000000;
std::vector<double> ua0, ua1, ub0, ub1, vola, volb;
}


void x_bar_makebar();
void x_bar_barcalc();


void x_bar(int, char**)
{
   char string[MAX_NCHAR];
   int exist = false;
   auto out = stdout;


   initial();


   int mode = 0;
   auto invalid_mode = [](int m) { return m != 1 and m != 2; };
   const char* mode_string1 = R"(
 The Tinker Thermodynamic Perturbation Utility Can :

    (1) Create BAR File with Perturbed Potential Energies
    (2) Compute Thermodynamic Values from Tinker BAR File
)";
   nextarg(string, exist);
   if (exist)
      read_string(mode, string);
   if (invalid_mode(mode))
      print(out, "%s", mode_string1);
   read_stream(mode,
               "\n"
               " Enter the Number of the Desired Choice :  ",
               0, invalid_mode);


   if (mode == 1)
      x_bar_makebar();
   if (mode == 2)
      x_bar_barcalc();
}


void x_bar_makebar()
{
   char string[MAX_NCHAR];
   int exist = false;


   auto invalid_temperature = [](double t) { return t <= 0; };
   auto invalid_recompute_answer = [](char c) {
      return c != 'y' and c != 'Y' and c != 'n' and c != 'N';
   };


   char answer;
   int recompute;
   int iarc, ibar;
   char arcfile[MAX_NCHAR];


   int lenga, lengb;
   int ltitlea, ltitleb;
   int nkey0, nkey1;
   int nfrma, nfrmb;
   double tempa, tempb;
   char filea[MAX_NCHAR], fileb[MAX_NCHAR];
   char titlea[MAX_NCHAR], titleb[MAX_NCHAR];
   std::vector<std::array<char, MAX_NCHAR>> keys0, keys1;


   // get trajectory A and setup mechanics calculation
   TINKER_RT(getarc)(&iarc);
   t_close(iarc);
   mechanic();


   // store the filename for trajectory A
   std::memcpy(filea, files::filename, MAX_NCHAR);
   lenga = files::leng;
   std::memcpy(titlea, titles::title, MAX_NCHAR);
   ltitlea = titles::ltitle;


   // store the keyword values for state 0
   nkey0 = keys::nkey;
   keys0.resize(nkey0);
   for (int i = 0; i < nkey0; ++i)
      std::memcpy(keys0[i].data(), keys::keyline[i], MAX_NCHAR);


   // find the original temperature value for trajectory A
   tempa = -1;
   nextarg(string, exist);
   if (exist)
      read_string(tempa, string);
   read_stream(tempa,
               "\n"
               " Enter Trajectory A Temperature in Degrees K [298] :  ",
               298.0, invalid_temperature);


   // get trajectory B and setup mechanics calculation
   TINKER_RT(getarc)(&iarc);
   t_close(iarc);
   mechanic();
   inform::silent = true;


   // store the filename for trajectory B
   std::memcpy(fileb, files::filename, MAX_NCHAR);
   lengb = files::leng;
   std::memcpy(titleb, titles::title, MAX_NCHAR);
   ltitleb = titles::ltitle;


   // store the keyword values for state 1
   nkey1 = keys::nkey;
   keys1.resize(nkey1);
   for (int i = 0; i < nkey1; ++i)
      std::memcpy(keys1[i].data(), keys::keyline[i], MAX_NCHAR);


   // find the original temperature value for trajectory B
   tempb = -1;
   nextarg(string, exist);
   if (exist)
      read_string(tempb, string);
   read_stream(tempb,
               "\n"
               " Enter Trajectory B Temperature in Degrees K [298] :  ",
               298.0, invalid_temperature);


   // decide whether to use energies from trajectory log files
   recompute = true;
   answer = ' ';
   nextarg(string, exist);
   if (exist)
      read_string(answer, string);
   read_stream(answer,
               "\n"
               " Obtain Energies from Trajectory Logs if Available [N] :  ",
               'N', invalid_recompute_answer);
   if (answer == 'y' or answer == 'Y')
      recompute = false;


   // reopen the file corresponding to trajectory A
   std::memcpy(arcfile, filea, MAX_NCHAR);
   t_suffix(arcfile, "arc", "old");
   std::string a_arc = arcfile;


   // reopen the file corresponding to trajectory B
   std::memcpy(arcfile, fileb, MAX_NCHAR);
   t_suffix(arcfile, "arc", "old");
   std::string b_arc = arcfile;


   // perform dynamic allocation of some local arrays
   ua0.resize(maxframe);
   ua1.resize(maxframe);
   ub0.resize(maxframe);
   ub1.resize(maxframe);
   vola.resize(maxframe);
   volb.resize(maxframe);


   //====================================================================//


   // check for log with energies of trajectory A in state 0


   // reset trajectory A using the parameters for state 0


   // find potential energies for trajectory A in state 0


   // reset trajectory A using the parameters for state 1


   // find potential energies for trajectory A in state 1


   // save potential energies and volumes for trajectory A


   //====================================================================//


   // check for log with energies of trajectory B in state 1


   // reset trajectory B using the parameters for state 1


   // find potential energies for trajectory B in state 1


   // reset trajectory B using the parameters for state 0


   // find potential energies for trajectory B in state 0


   // save potential energies and volumes for trajectory B


   //====================================================================//


   // perform deallocation of some local arrays
   ua0.clear();
   ua1.clear();
   ub0.clear();
   ub1.clear();
   vola.clear();
   volb.clear();
}


void x_bar_barcalc() {}
}
