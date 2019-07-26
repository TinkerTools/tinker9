#include "files.h"
#include "gpu/decl_basic.h"
#include "gpu/decl_tinker_io.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

using namespace TINKER_NAMESPACE;
using namespace test;

// static const int ans[5][216][70];
#include "test_nblist.hh"

static int usage_ =
    gpu::use_xyz | gpu::use_mass | gpu::use_traj | gpu::use_energy;

// if all integers in ref can be found in array
static bool find_match(const int* array, int na, int iframe, int iatom);

TEST_CASE("NBList-ArBox", "[ff][nblist][arbox]") {
  const char* k = "test_arbox.key";
  const char* x1 = "test_arbox.arc";
  const char* p = "amoeba09.prm";

  std::string k0 = arbox_key;
  file fke(k, k0);

  file fx1(x1, arbox_arc);
  file fpr(p, amoeba09_prm);

  const char* argv[] = {"dummy", x1};
  int argc = 2;
  test_begin_1_xyz(argc, argv);

  gpu::use_data = usage_;
  gpu::trajn = 5;
  tinker_gpu_data_create();

  gpu::copyin_tinker_arc(x1, 1, 5, 1);

  const int maxnlst = gpu::vlist_obj_.maxnlst;
  std::vector<int> nlst;
  nlst.resize(gpu::n);
  std::vector<int> lst;
  lst.resize(gpu::n * maxnlst);

  for (int ifr = 0;;) {
    gpu::copyout_array(nlst.data(), gpu::vlist_obj_.nlst, gpu::n);
    gpu::copyout_array(lst.data(), gpu::vlist_obj_.lst, gpu::n * maxnlst);

    for (int iatom = 0; iatom < gpu::n; ++iatom)
      REQUIRE(find_match(&lst[iatom * maxnlst], nlst[iatom], ifr, iatom));

    ++ifr;
    // update nblist
    if (ifr < gpu::trajn) {
      gpu::goto_frame0(ifr);
      gpu::evdw_reduce_xyz();
      gpu::nblist_data(rc_t::evolve);
    } else
      break;
  }

  tinker_gpu_data_destroy();
  test_end();
}

#include <set>
static bool find_match(const int* array, int na, int iframe, int iatom) {
  std::set<int> array_set(array, array + na);
  std::vector<int> missing;

  for (int i = 0; i < nmax; ++i) {
    int r = ans[iframe][iatom][i];
    if (r >= 0) {
      if (!array_set.count(r)) {
        missing.push_back(r);
      }
    }
  }

  if (missing.size() == 0) {
    return true;
  } else {
    std::string line = format(
        "The following neighbors of atom {:3d} in frame {:d} are missing:",
        iatom + 1, iframe + 1);
    for (auto i : missing) {
      line += format(" {:3d}", i + 1);
    }
    line += "\n";
    print(stdout, line);
    return false;
  }
}
