#include "test_rt.h"
#include "tinker_rt.h"
#include <fstream>
#include <tinker/detail/bath.hh>


namespace tinker {
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


void TestFile::keep()
{
   good_ = false;
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
#if TINKER_REAL_SIZE == 4
   (void)eps_double;
   return eps_single;
#elif TINKER_REAL_SIZE == 8
   (void)eps_single;
   return eps_double;
#else
   static_assert(false, "");
#endif
}


void test_begin_with_args(int argc, const char** argv)
{
   fortran_runtime_initialize(argc, const_cast<char**>(argv));

   initial();
   TINKER_RT(command)();
   TINKER_RT(getxyz)();
   mechanic();
   mechanic2();
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
}
