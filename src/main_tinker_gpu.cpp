#include "rc_man.h"
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

using namespace TINKER_NAMESPACE;

static const char* main_name = "tinker.gpu";
static const std::string analyze_name = "analyze";
static const std::string dynamic_name = "dynamic";
static const std::string helper_name = "help";
static const std::string info_name = "info";
static const std::string testgrad_name = "testgrad";
static const std::map<std::string, std::function<void(int, char**)>>&
launcher();

int main(int argc, char** argv)
{
   if (argc < 2) {
      goto help_message;
   }
   argc--;
   argv++;

   if (launcher().find(argv[0]) == launcher().end()) {
      goto help_message;
   } else {
      fortran_runtime_initialize(argc, argv);
      launcher().at(argv[0])(argc, argv);
      fortran_runtime_finish();
      return 0;
   }

help_message:
   launcher().at(helper_name)(argc, argv);
   std::exit(1);
}

TINKER_NAMESPACE_BEGIN
extern void x_analyze(int, char**);
extern void x_dynamic(int, char**);
extern void x_info(int, char**);
extern void x_testgrad(int, char**);
TINKER_NAMESPACE_END

static void x_help(int, char**)
{
   std::vector<std::string> keys;
   for (const auto& x : launcher()) {
      if (x.first != helper_name) {
         keys.push_back(x.first);
      }
   }
   std::sort(keys.begin(), keys.end());
   keys.push_back(helper_name);

   fprintf(stdout, " NAME\n");
   fprintf(stdout, "       %s\n\n", main_name);

   fprintf(stdout, " SYNOPSIS\n");
   fprintf(stdout, "       %s PROGRAM [ args... ]\n\n", main_name);

   fprintf(stdout, " PROGRAMS AVAILABLE\n");
   for (const auto& k : keys) {
      fprintf(stdout, "       %s\n", k.c_str());
   }
}

static const std::map<std::string, std::function<void(int, char**)>>& launcher()
{
   static std::map<std::string, std::function<void(int, char**)>> x = {
      {analyze_name, x_analyze},   {dynamic_name, x_dynamic},
      {helper_name, x_help},       {info_name, x_info},
      {testgrad_name, x_testgrad},
   };
   return x;
}
