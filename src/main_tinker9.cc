#include "main_tinker9.h"
#include "tinker_rt.h"
#include "tool/io_print.h"
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace tinker {
static const char* main_name = "tinker9";
static const std::string analyze_name = "analyze";
static const std::string bar_name = "bar";
static const std::string dynamic_name = "dynamic";
static const std::string helper_name = "help";
static const std::string info_name = "info";
static const std::string minimize_name = "minimize";
static const std::string testgrad_name = "testgrad";
static const std::map<std::string, std::function<void(int, char**)>>& launcher()
{
   static std::map<std::string, std::function<void(int, char**)>> x = {
      {analyze_name, x_analyze},
      {bar_name, x_bar},
      {dynamic_name, x_dynamic},
      {helper_name, x_help},
      {info_name, x_info},
      {minimize_name, x_minimize},
      {testgrad_name, x_testgrad},
   };
   return x;
}

static std::string getArg0(const char* a0)
{
   std::string s0 = a0;
   auto pos = s0.find_last_of('/') + 1;
   s0 = s0.substr(pos);
   const auto& l = launcher();
   std::string k;
   for (const auto& x : l) {
      k = x.first;
      if (k != helper_name)
         if (k + "9" == s0)
            return k;
   }
   return main_name;
}
}

int main(int argc, char** argv)
{
   using namespace tinker;

   std::string arg0;
   arg0 = getArg0(argv[0]);
   if (arg0 == main_name) {
      if (argc < 2)
         goto help_message;

      argc--;
      argv++;
      arg0 = argv[0];
   }

   if (launcher().find(arg0) == launcher().end()) {
      goto help_message;
   } else {
      tinkerFortranRuntimeBegin(argc, argv);
      try {
         launcher().at(arg0)(argc, argv);
      } catch (const std::exception& err) {
         print(stdout, " Terminating with uncaught exception :  %s\n",
               err.what());
      }
      tinkerFortranRuntimeEnd();
      return 0;
   }

help_message:
   launcher().at(helper_name)(argc, argv);
   std::exit(1);
}

namespace tinker {
void x_help(int, char**)
{
   std::vector<std::string> keys;
   for (const auto& x : launcher()) {
      if (x.first != helper_name) {
         keys.push_back(x.first);
      }
   }
   std::sort(keys.begin(), keys.end());
   keys.push_back(helper_name);

   print(stdout, " NAME\n");
   print(stdout, "       %s\n\n", main_name);

   print(stdout, " SYNOPSIS\n");
   print(stdout, "       %s PROGRAM [ args... ]\n\n", main_name);

   print(stdout, " PROGRAMS AVAILABLE\n");
   for (const auto& k : keys) {
      print(stdout, "       %s\n", k.c_str());
   }
}
}
