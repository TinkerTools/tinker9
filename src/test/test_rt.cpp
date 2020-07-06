#include "test_rt.h"
#include "tinker_rt.h"
#include "tool/error.h"
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


TestReference::TestReference(std::string fname)
{
   std::ifstream fr(fname);
   if (!fr)
      TINKER_THROW(format("TestReference cannot open file %s", fname));


   std::string l;
   while (fr) {
      std::getline(fr, l);
      Text::upcase(l);
      size_t end = std::string::npos;


      if (l.find("TOTAL POTENTIAL ENERGY :") != end) {
         //  Total Potential Energy :               -863.8791 Kcal/mole
         auto vs = Text::split(l);
         energy = std::stod(vs.end()[-2]);
      } else if (l.find("ENERGY COMPONENT BREAKDOWN :") != end) {
         std::getline(fr, l);
         auto vs = Text::split(l);
         energy = std::stod(vs.end()[-2]);
         count = std::stoi(vs.end()[-1]);
      } else if (l.find("INTERNAL VIRIAL TENSOR :") != end) {
         auto vs = Text::split(l);
         virial[0][0] = std::stod(vs.end()[-3]);
         virial[0][1] = std::stod(vs.end()[-2]);
         virial[0][2] = std::stod(vs.end()[-1]);
         std::getline(fr, l);
         vs = Text::split(l);
         virial[1][0] = std::stod(vs.end()[-3]);
         virial[1][1] = std::stod(vs.end()[-2]);
         virial[1][2] = std::stod(vs.end()[-1]);
         std::getline(fr, l);
         vs = Text::split(l);
         virial[2][0] = std::stod(vs.end()[-3]);
         virial[2][1] = std::stod(vs.end()[-2]);
         virial[2][2] = std::stod(vs.end()[-1]);
      } else if (l.find("ANLYT ") != end) {
         auto vs = Text::split(l);
         double g1 = std::stod(vs[2]);
         double g2 = std::stod(vs[3]);
         double g3 = std::stod(vs[4]);
         gradient.push_back(g1);
         gradient.push_back(g2);
         gradient.push_back(g3);
      }
   }
}


int TestReference::get_count() const
{
   return this->count;
}


double TestReference::get_energy() const
{
   return this->energy;
}


double (*TestReference::get_virial())[3]
{
   return this->virial;
}


double (*TestReference::get_gradient())[3]
{
   return reinterpret_cast<double(*)[3]>(this->gradient.data());
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
