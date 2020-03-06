#include "files.h"
#include "test.h"
#include "test_rt.h"


using namespace TINKER_NAMESPACE;


TEST_CASE("Vdw14-Trpcage", "[ff][evdw][lj][trpcage]")
{
   rc_flag = calc::xyz | calc::vmask;


   const char* kname = "test_vdw14.key";
   std::string k0 = trpcage_charmm19_key;
   k0 += "\nVDWTERM ONLY\n";
   const char* xname = "test_vdw14.xyz";
   const char* x0 = trpcage_charmm19_xyz;
   const char* pname = "charmm19.prm";
   const char* p0 = commit_11e84c69::charmm19_prm;


   const char* argv[] = {"dummy", x0};
   int argc = 2;


   SECTION("  - elj -- no pbc, no cutoff") {}


   SECTION("  - elj -- pbc, cutoff")
   {
      std::string k1 = k0 +
         "\nNEIGHBOR-LIST    0.5"
         "\nLIST-BUFFER      9.0"
         "\nA-AXIS            30"
         "\nB-AXIS            25"
         "\nC-AXIS            20"
         "\n";
   }
}
