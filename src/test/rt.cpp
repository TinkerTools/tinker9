#include "test/rt.h"
#include "util_rc_man.h"
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

file_expected::file_expected(const std::string& name) : name_(name) {}

file_expected::~file_expected() {
  std::ifstream chk(name_);
  if (chk) {
    std::remove(name_.c_str());
  }
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
}

void test_end() {
  TINKER_RT(final)();
  fortran_runtime_finish();
}

void test_mdinit(double t, double atm) {
  if (t > 0) {
    bath::kelvin = t;
    bath::isothermal = _true_;
  } else
    bath::isothermal = _false_;

  if (atm > 0) {
    bath::atmsph = atm;
    bath::isobaric = _true_;
  } else
    bath::isobaric = _false_;

  TINKER_RT(mdinit)();
}
}
TINKER_NAMESPACE_END
