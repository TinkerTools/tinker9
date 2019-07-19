#include "tinker_gpu.h"
#include <algorithm>
#include <functional>
#include <map>

m_tinker_using_namespace;

static const char* main_name = "tinker.gpu";
static const std::string dynamic_name = "dynamic";
static const std::string helper_name = "help";
static const std::map<std::string, std::function<void(int, char**)>>&
launcher();

int main(int argc, char** argv) {
  if (argc < 2) {
    goto help_message;
  }
  argc--;
  argv++;

  if (launcher().find(argv[0]) == launcher().end()) {
    fprintf(stdout, " Unrecognized program: %s\n\n", argv[0]);
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
extern void dynamic_x(int, char**);
TINKER_NAMESPACE_END

static void help_x(int, char**) {
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

static const std::map<std::string, std::function<void(int, char**)>>&
launcher() {
  static std::map<std::string, std::function<void(int, char**)>> x = {
      {dynamic_name, dynamic_x},
      {helper_name, help_x},
  };
  return x;
}
