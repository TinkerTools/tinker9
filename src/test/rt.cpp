#include "test/rt.h"
#include "util/fort_rt.h"
#include <cstdio>
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>
#include <fstream>

TINKER_NAMESPACE_BEGIN
namespace test {
file::file(const std::string& name, const std::string& content)
    : good_(false), name_(name) {
  std::ofstream fout(name);
  good_ = fout.is_open();
  if (good_) {
    fout << content;
  }
}

file::~file() {
  if (good_)
    std::remove(name_.c_str());
}

double test_get_eps2(double eps_single, double eps_double) {
#if defined(TINKER_GPU_SINGLE)
  return eps_single;
#elif defined(TINKER_GPU_DOUBLE)
  return eps_double;
#else
  static_assert(false, "");
#endif
}

void test_begin_1_xyz(int argc, const char** argv) {
  fortran_runtime_initialize(argc, const_cast<char**>(argv));

  TINKER_RT(initial)();
  TINKER_RT(command)();
  TINKER_RT(getxyz)();
  TINKER_RT(mechanic)();

  bath::isothermal = _true_;
  bath::kelvin = 298;
  bath::isobaric = _true_;
  bath::atmsph = 1;

  TINKER_RT(mdinit)();
}

void test_end() {
  TINKER_RT(final)();
  fortran_runtime_finish();
}
}
TINKER_NAMESPACE_END
