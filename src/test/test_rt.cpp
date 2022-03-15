#include "test_rt.h"
#include "tinker_rt.h"
#include "tool/error.h"
#include <fstream>
#include <tinker/detail/bath.hh>

namespace tinker {
TestFile::TestFile(const std::string& file, std::string dst, std::string extra)
{
   if (dst == "") {
      auto pos = file.find_last_of('/');
      name = file.substr(pos + 1);
      if (name == "")
         return;
   } else {
      name = dst;
   }

   if (file != "") {
      std::ifstream fsrc(file, std::ios::binary);
      std::ofstream fdst(name, std::ios::binary);
      good = fdst.is_open();
      if (good) {
         fdst << fsrc.rdbuf();
         if (extra != "") {
            fdst << extra;
         }
      }
   } else {
      std::ofstream fout(name);
      good = fout.is_open();
      if (good) {
         fout << extra;
      }
   }
}

TestFile::~TestFile()
{
   if (good)
      std::remove(name.c_str());
}

void TestFile::__keep()
{
   good = false;
}

TestFileExpected::TestFileExpected(const std::string& name)
   : m_name(name)
{}

TestFileExpected::~TestFileExpected()
{
   std::ifstream chk(m_name);
   if (chk) {
      std::remove(m_name.c_str());
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
   tinkerFortranRuntimeBegin(argc, const_cast<char**>(argv));

   initial();
   tinker_f_command();
   tinker_f_getxyz();
   tinker_f_mechanic();
   mechanic2();
}

void test_end()
{
   tinker_f_final();
   tinkerFortranRuntimeEnd();
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

   double dt = 0.001;
   tinker_f_mdinit(&dt);
}
}
