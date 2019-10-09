#include "test_rt.h"
#include "rc_man.h"
#include "tinker_rt.h"
#include <ext/tinker/detail/bath.hh>
#include <fstream>

TINKER_NAMESPACE_BEGIN
TestFile::TestFile(const std::string& name, const std::string& content)
   : good_(false)
   , name_(name)
{
   std::ofstream fout(name);
   good_ = fout.is_open();
   if (good_) {
      fout << content;
   }
}

TestFile::~TestFile()
{
   if (good_)
      std::remove(name_.c_str());
}

TestFileExpected::TestFileExpected(const std::string& name)
   : name_(name)
{}

TestFileExpected::~TestFileExpected()
{
   std::ifstream chk(name_);
   if (chk) {
      std::remove(name_.c_str());
   }
}

double test_get_eps(double eps_single, double eps_double)
{
#if TINKER_SINGLE_PRECISION
   return eps_single;
#elif TINKER_DOUBLE_PRECISION
   return eps_double;
#else
   static_assert(false, "");
#endif
}

void test_begin_with_args(int argc, const char** argv)
{
   fortran_runtime_initialize(argc, const_cast<char**>(argv));

   TINKER_RT(initial)();
   TINKER_RT(command)();
   TINKER_RT(getxyz)();
   TINKER_RT(mechanic)();
}

void test_end()
{
   TINKER_RT(final)();
   fortran_runtime_finish();
}

void test_mdinit(double t, double atm)
{
   if (t > 0) {
      bath::kelvin = t;
      bath::isothermal = true;
   } else
      bath::isothermal = false;

   if (atm > 0) {
      bath::atmsph = atm;
      bath::isobaric = true;
   } else
      bath::isobaric = false;

   TINKER_RT(mdinit)();
}
TINKER_NAMESPACE_END
