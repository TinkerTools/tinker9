#include "tool/error.h"
#include <tinker/detail/bath.hh>
#include <tinker/routines.h>

#include "testrt.h"
#include "tinker9.h"

#include <fstream>
#include <map>
#include <tuple>

namespace tinker {
TestFile::TestFile(std::string file, std::string dst, std::string extra)
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

TestRemoveFileOnExit::TestRemoveFileOnExit(std::string name)
   : m_name(name)
{}

TestRemoveFileOnExit::~TestRemoveFileOnExit()
{
   std::ifstream chk(m_name);
   if (chk) {
      std::remove(m_name.c_str());
   }
}
}

namespace tinker {
class TestReference::Impl
{
public:
   std::vector<double> gradient;
   std::map<std::string, std::tuple<double, int>> engcnt;
   double virial[3][3];
   double energy;
   int count;
};

TestReference::~TestReference()
{
   delete pimpl;
}

TestReference::TestReference(std::string fname)
   : pimpl(new TestReference::Impl)
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
         pimpl->energy = std::stod(vs.end()[-2]);
      } else if (l.find("ENERGY COMPONENT BREAKDOWN :") != end) {
         std::getline(fr, l);
         auto vs = Text::split(l);
         pimpl->energy = std::stod(vs.end()[-2]);
         pimpl->count = std::stoi(vs.end()[-1]);
      } else if (l.find("INTERNAL VIRIAL TENSOR :") != end) {
         auto vs = Text::split(l);
         pimpl->virial[0][0] = std::stod(vs.end()[-3]);
         pimpl->virial[0][1] = std::stod(vs.end()[-2]);
         pimpl->virial[0][2] = std::stod(vs.end()[-1]);
         std::getline(fr, l);
         vs = Text::split(l);
         pimpl->virial[1][0] = std::stod(vs.end()[-3]);
         pimpl->virial[1][1] = std::stod(vs.end()[-2]);
         pimpl->virial[1][2] = std::stod(vs.end()[-1]);
         std::getline(fr, l);
         vs = Text::split(l);
         pimpl->virial[2][0] = std::stod(vs.end()[-3]);
         pimpl->virial[2][1] = std::stod(vs.end()[-2]);
         pimpl->virial[2][2] = std::stod(vs.end()[-1]);
      } else if (l.find("ANLYT ") != end) {
         auto vs = Text::split(l);
         double g1 = std::stod(vs[2]);
         double g2 = std::stod(vs[3]);
         double g3 = std::stod(vs[4]);
         pimpl->gradient.push_back(g1);
         pimpl->gradient.push_back(g2);
         pimpl->gradient.push_back(g3);
      } else if (l.find("ENGCNT ") != end) {
         auto vs = Text::split(l);
         std::string name = vs[1];
         for (size_t i = 2; i + 2 < vs.size(); ++i) {
            name += " ";
            name += vs[i];
         }
         double eng = std::stod(vs.end()[-2]);
         int cnt = std::stoi(vs.end()[-1]);
         pimpl->engcnt[name] = std::make_tuple(eng, cnt);
      }
   }
}

int TestReference::getCount() const
{
   return pimpl->count;
}

double TestReference::getEnergy() const
{
   return pimpl->energy;
}

void TestReference::getEnergyCountByName(std::string name, double& energy, int& count)
{
   Text::upcase(name);
   auto& it = pimpl->engcnt.at(name);
   energy = std::get<0>(it);
   count = std::get<1>(it);
}

const double (*TestReference::getVirial() const)[3]
{
   return pimpl->virial;
}

const double (*TestReference::getGradient() const)[3]
{
   return reinterpret_cast<const double(*)[3]>(pimpl->gradient.data());
}

double testGetEps(double eps_single, double eps_double)
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

void testBeginWithArgs(int argc, const char** argv)
{
   tinkerFortranRuntimeBegin(argc, const_cast<char**>(argv));

   initial();
   tinker_f_command();
   tinker_f_getxyz();
   tinker_f_mechanic();
   mechanic2();
}

void testEnd()
{
   tinker_f_final();
   tinkerFortranRuntimeEnd();
}

void testMdInit(double t, double atm)
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
