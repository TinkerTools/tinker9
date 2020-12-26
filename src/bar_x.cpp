#include "energy.h"
#include "nblist.h"
#include "tinker_rt.h"
#include <array>
#include <fstream>
#include <tinker/detail/files.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/keys.hh>
#include <tinker/detail/titles.hh>


namespace tinker {
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


   TINKER_RT(final)();
}


void x_bar_makebar()
{
   rc_flag = calc::xyz | calc::mass | calc::energy;


   auto out = stdout;
   auto invalid_temperature = [](double t) { return t <= 0; };
   auto invalid_recompute_answer = [](char c) {
      return c != 'y' and c != 'Y' and c != 'n' and c != 'N';
   };


   int lenga, lengb;
   int ltitlea, ltitleb;
   int nkey0, nkey1;
   double tempa, tempb;
   char filea[MAX_NCHAR], fileb[MAX_NCHAR];
   char titlea[MAX_NCHAR], titleb[MAX_NCHAR];
   std::vector<std::array<char, MAX_NCHAR>> keys0, keys1;
   char string[MAX_NCHAR];
   std::string str;


   // get trajectory A and setup mechanics calculation
   int iarc;
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
   int exist = false;
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
   int recompute = true;
   char answer = ' ';
   nextarg(string, exist);
   if (exist)
      read_string(answer, string);
   read_stream(answer,
               "\n"
               " Obtain Energies from Trajectory Logs if Available [N] :  ",
               'N', invalid_recompute_answer);
   if (answer == 'y' or answer == 'Y')
      recompute = false;


   //====================================================================//


   int done;
   const char* log_fmt = " Current Potential %lf\n";
   const char* process_fmt = "       Completed%8d Coordinate Frames\n";
   std::memcpy(string, filea, lenga);
   t_suffix(string, "bar", "new");
   std::string barfile = string;
   std::vector<double> ua0, ua1, ub0, ub1, vola, volb;


   // check for log with energies of trajectory A in state 0
   if (not recompute) {
      std::memcpy(string, filea, lenga);
      t_suffix(string, "log", "old");
      str = string;
      std::ifstream a_log(str);
      while (std::getline(a_log, str)) {
         double val;
         int count = std::sscanf(str.c_str(), log_fmt, &val);
         if (count == 1)
            ua0.push_back(val);
      }
      a_log.close();
   }


   // reopen the file corresponding to trajectory A
   print(out,
         "\n"
         " Initial Processing for Trajectory A :\n"
         "\n");


   std::memcpy(string, filea, MAX_NCHAR);
   t_suffix(string, "arc", "old");
   iarc = t_freeunit();
   t_open(iarc, string, "old");
   str = string;
   std::ifstream a_arc(str);


   if (ua0.size() == 0) {
      // reset trajectory A using the parameters for state 0
      rewind_stream(a_arc);
      t_rewind(iarc);
      TINKER_RT(readxyz)(&iarc);
      keys::nkey = nkey0;
      for (int i = 0; i < keys::nkey; ++i)
         std::memcpy(keys::keyline[i], keys0[i].data(), MAX_NCHAR);
      mechanic();


      // find potential energies for trajectory A in state 0
      initialize();
      done = false;
      do {
         read_frame_copyin_to_xyz(a_arc, done);
         refresh_neighbors();
         energy(calc::energy);
         ua0.push_back(esum);
         int i = ua0.size();
         if (i % 100 == 0) {
            print(out, process_fmt, i);
            std::fflush(out);
         }
      } while (not done);
      finish();
   } else {
      int ii = ua0.size();
      for (int i = 1; i <= ii; ++i) {
         if (i % 100 == 0)
            print(out, process_fmt, i);
      }
      std::fflush(out);
   }


   // reset trajectory A using the parameters for state 1
   rewind_stream(a_arc);
   t_rewind(iarc);
   TINKER_RT(readxyz)(&iarc);
   keys::nkey = nkey1;
   for (int i = 0; i < keys::nkey; ++i)
      std::memcpy(keys::keyline[i], keys1[i].data(), MAX_NCHAR);
   mechanic();


   // find potential energies for trajectory A in state 1
   if (inform::verbose)
      print(out,
            "\n"
            " Potential Energy Values for Trajectory A :"
            "\n\n"
            "       Frame         State 0         State 1            Delta\n"
            "\n");
   initialize();
   done = false;
   do {
      read_frame_copyin_to_xyz(a_arc, done);
      refresh_neighbors();
      energy(calc::energy);
      double vol = volbox();
      ua1.push_back(esum);
      vola.push_back(vol);
      if (inform::verbose) {
         int i = ua1.size() - 1;
         print(out, "%11d  %16.4lf%16.4lf%16.4lf\n", i + 1, ua0[i], ua1[i],
               ua1[i] - ua0[i]);
      }
   } while (not done);
   finish();


   // save potential energies and volumes for trajectory A
   t_close(iarc);
   FILE* ibar = std::fopen(barfile.c_str(), "w");
   int nfrma = std::min(ua0.size(), ua1.size());
   str = fstr_view(titlea)(1, ltitlea).trim();
   print(ibar, "%8d%10.2lf  %s\n", nfrma, tempa, str);
   for (int i = 0; i < nfrma; ++i) {
      if (vola[i] == 0)
         print(ibar, "%8d  %18.4lf%18.4lf\n", i + 1, ua0[i], ua1[i]);
      else
         print(ibar, "%8d  %18.4lf%18.4lf%18.4lf\n", i + 1, ua0[i], ua1[i],
               vola[i]);
   }
   std::fflush(ibar);


   //====================================================================//


   // check for log with energies of trajectory B in state 1
   if (not recompute) {
      std::memcpy(string, fileb, lengb);
      t_suffix(string, "log", "old");
      str = string;
      std::ifstream b_log(str);
      while (std::getline(b_log, str)) {
         double val;
         int count = std::sscanf(str.c_str(), log_fmt, &val);
         if (count == 1)
            ub1.push_back(val);
      }
      b_log.close();
   }


   // reopen the file corresponding to trajectory B
   print(out,
         "\n"
         " Initial Processing for Trajectory B :\n"
         "\n");


   std::memcpy(string, fileb, MAX_NCHAR);
   t_suffix(string, "arc", "old");
   iarc = t_freeunit();
   t_open(iarc, string, "old");
   str = string;
   std::ifstream b_arc(str);


   if (ub1.size() == 0) {
      // reset trajectory B using the parameters for state 1
      rewind_stream(b_arc);
      t_rewind(iarc);
      TINKER_RT(readxyz)(&iarc);
      keys::nkey = nkey1;
      for (int i = 0; i < keys::nkey; ++i)
         std::memcpy(keys::keyline[i], keys1[i].data(), MAX_NCHAR);
      mechanic();


      // find potential energies for trajectory B in state 1
      initialize();
      done = false;
      do {
         read_frame_copyin_to_xyz(b_arc, done);
         refresh_neighbors();
         energy(calc::energy);
         ub1.push_back(esum);
         int i = ub1.size();
         if (i % 100 == 0) {
            print(out, process_fmt, i);
            std::fflush(out);
         }
      } while (not done);
      finish();
   } else {
      int ii = ub1.size();
      for (int i = 1; i <= ii; ++i) {
         if (i % 100 == 0)
            print(out, process_fmt, i);
      }
      std::fflush(out);
   }


   // reset trajectory B using the parameters for state 0
   rewind_stream(b_arc);
   t_rewind(iarc);
   TINKER_RT(readxyz)(&iarc);
   keys::nkey = nkey0;
   for (int i = 0; i < keys::nkey; ++i)
      std::memcpy(keys::keyline[i], keys0[i].data(), MAX_NCHAR);
   mechanic();


   // find potential energies for trajectory B in state 0
   if (inform::verbose)
      print(out,
            "\n"
            " Potential Energy Values for Trajectory B :"
            "\n\n"
            "       Frame         State 0         State 1            Delta\n"
            "\n");
   initialize();
   done = false;
   do {
      read_frame_copyin_to_xyz(b_arc, done);
      refresh_neighbors();
      energy(calc::energy);
      double vol = volbox();
      ub0.push_back(esum);
      volb.push_back(vol);
      if (inform::verbose) {
         int i = ub0.size() - 1;
         print(out, "%11d  %16.4lf%16.4lf%16.4lf\n", i + 1, ub0[i], ub1[i],
               ub0[i] - ub1[i]);
      }
   } while (not done);
   finish();


   // save potential energies and volumes for trajectory B
   t_close(iarc);
   int nfrmb = std::min(ub0.size(), ub1.size());
   str = fstr_view(titleb)(1, ltitleb).trim();
   print(ibar, "%8d%10.2lf  %s\n", nfrmb, tempb, str);
   for (int i = 0; i < nfrmb; ++i) {
      if (volb[i] == 0)
         print(ibar, "%8d  %18.4lf%18.4lf\n", i + 1, ub0[i], ub1[i]);
      else
         print(ibar, "%8d  %18.4lf%18.4lf%18.4lf\n", i + 1, ub0[i], ub1[i],
               volb[i]);
   }
   std::fclose(ibar);
   print(out,
         "\n"
         " Potential Energy Values Written To :  %s\n",
         barfile);
}


void x_bar_barcalc()
{
   auto out = stdout;
   int exist = false;
   int query;
   char string[MAX_NCHAR];
   std::string str;
   auto invalid_barfile = [&](const std::string& s) {
      std::memcpy(string, s.data(), s.length());
      t_basefile(string);
      t_suffix(string, "bar", "old");
      std::ifstream fr(s);
      return not fr;
   };


   // ask the user for file with potential energies and volumes
   nextarg(string, exist);
   if (exist) {
      t_basefile(string);
      t_suffix(string, "bar", "old");
      str = string;
   }
   read_stream(str,
               "\n"
               " Enter Potential Energy BAR File Name :  ",
               std::string(""), invalid_barfile);
   std::string barfile = string;


   auto invalid_frames = [](const std::string& s) {
      std::istringstream iss(s);
      int start, stop, step;
      iss >> start >> stop >> step;
      return start <= 0 or stop <= start or step <= 0;
   };


   // set beginning and ending frame for trajectory A
   int starta = 0, stopa = 0, stepa = 0;
   query = true;
   if (query) {
      nextarg(string, exist);
      if (exist) {
         int count = std::sscanf(string, "%d", &starta);
         query = count == 1;
      }
   }
   if (query) {
      nextarg(string, exist);
      if (exist) {
         int count = std::sscanf(string, "%d", &stopa);
         query = count == 1;
      }
   }
   if (query) {
      nextarg(string, exist);
      if (exist) {
         int count = std::sscanf(string, "%d", &stepa);
         query = count == 1;
      }
   }
   str = std::to_string(starta) + " " + std::to_string(stopa) + " " +
      std::to_string(stepa);
   read_stream(str,
               "\n"
               " First & Last Frame and Step Increment for Trajectory A :  ",
               std::string(""), invalid_frames);
   {
      std::istringstream iss(str);
      iss >> starta >> stopa >> stepa;
   }


   // set beginning and ending frame for trajectory B
   int startb = 0, stopb = 0, stepb = 0;
   query = true;
   if (query) {
      nextarg(string, exist);
      if (exist) {
         int count = std::sscanf(string, "%d", &startb);
         query = count == 1;
      }
   }
   if (query) {
      nextarg(string, exist);
      if (exist) {
         int count = std::sscanf(string, "%d", &stopb);
         query = count == 1;
      }
   }
   if (query) {
      nextarg(string, exist);
      if (exist) {
         int count = std::sscanf(string, "%d", &stepb);
         query = count == 1;
      }
   }
   str = std::to_string(startb) + " " + std::to_string(stopb) + " " +
      std::to_string(stepb);
   str = std::to_string(startb) + " " + std::to_string(stopb) + " " +
      std::to_string(stepb);
   read_stream(str,
               "\n"
               " First & Last Frame and Step Increment for Trajectory B :  ",
               std::string(""), invalid_frames);
   {
      std::istringstream iss(str);
      iss >> startb >> stopb >> stepb;
   }


   std::vector<double> ua0, ua1, ub0, ub1, vola, volb, vloga, vlogb;
   std::ifstream ibar(barfile);


   //====================================================================//


   // read potential energies and volumes for trajectory A
   int nfrma = 0, nfrma1;
   double tempa = 0;
   std::getline(ibar, str);
   std::sscanf(str.c_str(), "%d %lf\n", &nfrma, &tempa);
   nfrma1 = nfrma;
   str = str.substr(str.find_first_not_of(' '));
   str = str.substr(str.find_first_of(' '));
   str = str.substr(str.find_first_not_of(' '));
   str = str.substr(str.find_first_of(' '));
   str = str.substr(str.find_first_not_of(' '));
   std::string titlea = str;
   nfrma = std::min(nfrma, stopa);
   for (int i = 1; i <= nfrma1; ++i) {
      std::getline(ibar, str);
      if (((i - starta) % stepa) == 0 and i <= stopa) {
         double u0 = 0, u1 = 0, v = 0;
         std::sscanf(str.c_str(), "%*d %lf %lf %lf\n", &u0, &u1, &v);
         ua0.push_back(u0);
         ua1.push_back(u1);
         vola.push_back(v);
      }
   }
   nfrma = vola.size();


   // read potential energies and volumes for trajectory B
   int nfrmb = 0, nfrmb1;
   double tempb = 0;
   std::getline(ibar, str);
   std::sscanf(str.c_str(), "%d %lf\n", &nfrmb, &tempb);
   nfrmb1 = nfrmb;
   str = str.substr(str.find_first_not_of(' '));
   str = str.substr(str.find_first_of(' '));
   str = str.substr(str.find_first_not_of(' '));
   str = str.substr(str.find_first_of(' '));
   str = str.substr(str.find_first_not_of(' '));
   std::string titleb = str;
   nfrmb = std::min(nfrmb, stopb);
   for (int i = 1; i <= nfrmb1; ++i) {
      std::getline(ibar, str);
      if (((i - startb) % stepb) == 0 and i <= stopb) {
         double u0 = 0, u1 = 0, v = 0;
         std::sscanf(str.c_str(), "%*d %lf %lf %lf\n", &u0, &u1, &v);
         ub0.push_back(u0);
         ub1.push_back(u1);
         volb.push_back(v);
      }
   }
   nfrmb = volb.size();
   ibar.close();


   // provide info about trajectories and number of frames
   print(out,
         "\n"
         " Simulation Trajectory A and Thermodynamic State 0 :  \n"
         "\n"
         " %s\n",
         titlea);
   print(out,
         " Number of Frames :    %8d\n"
         " Temperature :       %10.2lf\n",
         nfrma, tempa);
   print(out,
         "\n"
         " Simulation Trajectory B and Thermodynamic State 1 :  \n"
         "\n"
         " %s\n",
         titleb);
   print(out,
         " Number of Frames :    %8d\n"
         " Temperature :       %10.2lf\n",
         nfrmb, tempb);
}
}
