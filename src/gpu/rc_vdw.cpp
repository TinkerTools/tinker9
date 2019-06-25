#include "gpu/decl_mdstate.h"
#include "gpu/decl_switch.h"
#include "gpu/e_vdw.h"
#include "gpu/rc.h"
#include "util/fort_str.h"
#include <map>

TINKER_NAMESPACE_BEGIN
namespace gpu {
evdw_t vdwtyp;

const char* vdwtyp_str(evdw_t typ) {
  if (typ == vdw_lj)
    return "LENNARD-JONES";
  else if (typ == vdw_buck)
    return "BUCKINGHAM";
  else if (typ == vdw_mm3hb)
    return "MM3-HBOND";
  else if (typ == vdw_hal)
    return "BUFFERED-14-7";
  else if (typ == vdw_gauss)
    return "GAUSSIAN";
  else {
    assert(false);
    return "";
  }
}

double vdw_switch_cut, vdw_switch_off;

int* ired;
real* kred;
real *xred, *yred, *zred;

int *jvdw, *njvdw;
real *radmin, *epsilon;

real* vlam;

real* ev;
int* nev;
real* vir_ev;
int use_evdw() { return potent::use_vdw; }

void get_evdw_type(evdw_t& typ) {
  fstr_view str = vdwpot::vdwtyp;
  if (str == "LENNARD-JONES")
    typ = vdw_lj;
  else if (str == "BUCKINGHAM")
    typ = vdw_buck;
  else if (str == "MM3-HBOND")
    typ = vdw_mm3hb;
  else if (str == "BUFFERED-14-7")
    typ = vdw_hal;
  else if (str == "GAUSSIAN")
    typ = vdw_gauss;
  else
    assert(false);
}

void evdw_data(rc_t rc) {
  if (!use_evdw())
    return;

  typedef int new_type;
  typedef int old_type;
  static std::map<old_type, new_type> jmap;
  static std::vector<new_type> jvec;
  static std::vector<new_type> jvdwbuf;
  static int jcount;

  if (rc & rc_dealloc) {
    // local static members
    jmap.clear();
    jvec.clear();
    jvdwbuf.clear();
    jcount = 0;

    check_cudart(cudaFree(ired));
    check_cudart(cudaFree(kred));
    check_cudart(cudaFree(xred));
    check_cudart(cudaFree(yred));
    check_cudart(cudaFree(zred));

    check_cudart(cudaFree(jvdw));
    check_cudart(cudaFree(njvdw));
    check_cudart(cudaFree(radmin));
    check_cudart(cudaFree(epsilon));

    check_cudart(cudaFree(vlam));

    check_cudart(cudaFree(ev));
    check_cudart(cudaFree(nev));
    check_cudart(cudaFree(vir_ev));
  }

  if (rc & rc_alloc) {
    const size_t rs = sizeof(real);
    size_t size;

    size = n * rs;
    check_cudart(cudaMalloc(&ired, n * sizeof(int)));
    check_cudart(cudaMalloc(&kred, size));
    check_cudart(cudaMalloc(&xred, size));
    check_cudart(cudaMalloc(&yred, size));
    check_cudart(cudaMalloc(&zred, size));

    check_cudart(cudaMalloc(&jvdw, n * sizeof(int)));
    check_cudart(cudaMalloc(&njvdw, sizeof(int)));

    jvdwbuf.resize(n);
    assert(jmap.size() == 0);
    assert(jvec.size() == 0);
    jcount = 0;
    for (int i = 0; i < n; ++i) {
      int jt = vdw::jvdw[i] - 1;
      auto iter = jmap.find(jt);
      if (iter == jmap.end()) {
        jvdwbuf[i] = jcount;
        jvec.push_back(jt);
        jmap[jt] = jcount;
        ++jcount;
      } else {
        jvdwbuf[i] = iter->second;
      }
    }
    size = jcount * jcount * rs;
    check_cudart(cudaMalloc(&radmin, size));
    check_cudart(cudaMalloc(&epsilon, size));

    size = n * rs;
    check_cudart(cudaMalloc(&vlam, size));

    check_cudart(cudaMalloc(&ev, rs));
    check_cudart(cudaMalloc(&nev, sizeof(int)));
    check_cudart(cudaMalloc(&vir_ev, rs * 9));
  }

  if (rc & rc_copyin) {
    get_evdw_type(vdwtyp);

    switch_cut_off(switch_vdw, vdw_switch_cut, vdw_switch_off);

    std::vector<int> iredbuf(n);
    std::vector<double> kredbuf(n);
    for (int i = 0; i < n; ++i) {
      int jt = vdw::ired[i] - 1;
      iredbuf[i] = jt;
      kredbuf[i] = vdw::kred[i];
    }
    copyin_array(ired, iredbuf.data(), n);
    copyin_array(kred, kredbuf.data(), n);

    copyin_array(jvdw, jvdwbuf.data(), n);
    copyin_array(njvdw, &jcount, 1);

    // see also kvdw.f
    std::vector<double> radvec, epsvec;
    for (int it_new = 0; it_new < jcount; ++it_new) {
      int it_old = jvec[it_new];
      int base = it_old * sizes::maxclass;
      for (int jt_new = 0; jt_new < jcount; ++jt_new) {
        int jt_old = jvec[jt_new];
        int offset = base + jt_old;
        radvec.push_back(vdw::radmin[offset]);
        epsvec.push_back(vdw::epsilon[offset]);
      }
    }
    copyin_array(radmin, radvec.data(), jcount * jcount);
    copyin_array(epsilon, epsvec.data(), jcount * jcount);

    std::vector<double> vlamvec(n);
    for (int i = 0; i < n; ++i) {
      if (mutant::mut[i]) {
        vlamvec[i] = mutant::vlambda;
      } else {
        vlamvec[i] = 1;
      }
    }
    copyin_array(vlam, vlamvec.data(), n);
  }
}

void evdw(int vers) {
  if (vdwtyp == vdw_lj)
    evdw_lj(vers);
  else if (vdwtyp == vdw_buck)
    evdw_buck(vers);
  else if (vdwtyp == vdw_mm3hb)
    evdw_mm3hb(vers);
  else if (vdwtyp == vdw_hal)
    evdw_hal(vers);
  else if (vdwtyp == vdw_gauss)
    evdw_gauss(vers);
}
}
TINKER_NAMESPACE_END
