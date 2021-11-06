#include "main_tinker9.h"
#include "tool/io_print.h"
#include <tinker/routines.h>

#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>


namespace tinker {
namespace {
const char* main_name = "tinker9";
const std::string analyze_name = "analyze";
const std::string bar_name = "bar";
const std::string dynamic_name = "dynamic";
const std::string helper_name = "help";
const std::string info_name = "info";
const std::string minimize_name = "minimize";
const std::string testgrad_name = "testgrad";
const std::map<std::string, std::function<void(int, char**)>>& launcher()
{
   static std::map<std::string, std::function<void(int, char**)>> x = {
      {analyze_name, x_analyze},   {bar_name, x_bar},
      {dynamic_name, x_dynamic},   {helper_name, x_help},
      {info_name, x_info},         {minimize_name, x_minimize},
      {testgrad_name, x_testgrad},
   };
   return x;
}
}


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


int main(int argc, char** argv)
{
   using namespace tinker;

   if (argc < 2) {
      goto help_message;
   }
   argc--;
   argv++;

   if (launcher().find(argv[0]) == launcher().end()) {
      goto help_message;
   } else {
      tinkerFortranRuntimeBegin(argc, argv);
      try {
         launcher().at(argv[0])(argc, argv);
      } catch (const std::exception& err) {
         fprintf(stdout, " Terminating with uncaught exception :  %s\n",
                 err.what());
      }
      tinkerFortranRuntimeEnd();
      return 0;
   }

help_message:
   launcher().at(helper_name)(argc, argv);
   std::exit(1);
}
