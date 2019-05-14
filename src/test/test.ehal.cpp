#include "files.h"
#include "test/rt.h"
#include "test/test.h"

m_tinker_using_namespace;
using namespace test;

TEST_CASE("EHal") {
  SECTION("1") {
    const char* xyz = "test_watersmall.xyz";
    const char* prm = "water03.prm";
    const char* key = "test_watersmall.key";

    const char* key_content = R"**(
parameters                  water03
verbose

a-axis                       18.643
integrator                    RESPA
thermostat                    BUSSI
barostat                 MONTECARLO

#neighbor-list
list-buffer                     0.1
vdw-cutoff                      9.0
ewald
ewald-cutoff                    7.0

polarization                 MUTUAL
polar-eps                   0.00001

vdwterm                        only
)**";

    file f_xyz(xyz, test::watersmall_xyz);
    file f_prm(prm, test::water03_prm);
    file f_key(key, key_content);

    std::string cmd = "analyze ";
    cmd += xyz;
    cmd += " e";
    int ret = std::system(cmd.c_str());
    REQUIRE(ret == 0);
  }
}
