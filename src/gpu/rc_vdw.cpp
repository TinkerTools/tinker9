#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"
#include "gpu/decl_switch.h"
#include "gpu/e_vdw.h"
#include "rc_cudart.h"
#include "util/fort_str.h"
#include <map>

TINKER_NAMESPACE_BEGIN
namespace gpu {
int vdwtyp = 0;
std::string vdwtyp_str;

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

void get_evdw_type(int& typ, std::string& typ_str) {
  fstr_view str = vdwpot::vdwtyp;
  typ_str = str.trim();
  if (str == "LENNARD-JONES")
    typ = evdw_lj;
  else if (str == "BUCKINGHAM")
    typ = evdw_buck;
  else if (str == "MM3-HBOND")
    typ = evdw_mm3hb;
  else if (str == "BUFFERED-14-7")
    typ = evdw_hal;
  else if (str == "GAUSSIAN")
    typ = evdw_gauss;
}

void evdw_data(int op) {
  if (!use_evdw())
    return;

  if (op & op_dealloc) {
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

  // TODO
  if (op & op_alloc) {
  }
  if (op & op_copyin) {
  }

  if (op == op_create) {
    get_evdw_type(vdwtyp, vdwtyp_str);

    const size_t rs = sizeof(real);
    size_t size;

    switch_cut_off(switch_vdw, vdw_switch_cut, vdw_switch_off);

    size = n * rs;
    check_cudart(cudaMalloc(&ired, n * sizeof(int)));
    check_cudart(cudaMalloc(&kred, size));
    check_cudart(cudaMalloc(&xred, size));
    check_cudart(cudaMalloc(&yred, size));
    check_cudart(cudaMalloc(&zred, size));
    std::vector<int> iredbuf(n);
    std::vector<double> kredbuf(n);
    for (int i = 0; i < n; ++i) {
      int jt = vdw::ired[i] - 1;
      iredbuf[i] = jt;
      kredbuf[i] = vdw::kred[i];
    }
    copyin_data(ired, iredbuf.data(), n);
    copyin_data(kred, kredbuf.data(), n);

    check_cudart(cudaMalloc(&jvdw, n * sizeof(int)));
    check_cudart(cudaMalloc(&njvdw, sizeof(int)));
    std::vector<int> jbuf(n);
    typedef int new_type;
    typedef int old_type;
    std::map<old_type, new_type> jmap;
    std::vector<new_type> jvec;
    int jcount = 0;
    for (int i = 0; i < n; ++i) {
      int jt = vdw::jvdw[i] - 1;
      auto iter = jmap.find(jt);
      if (iter == jmap.end()) {
        jbuf[i] = jcount;
        jvec.push_back(jt);
        jmap[jt] = jcount;
        ++jcount;
      } else {
        jbuf[i] = iter->second;
      }
    }
    copyin_data(jvdw, jbuf.data(), n);
    copyin_data(njvdw, &jcount, 1);
    size = jcount * jcount * rs;
    check_cudart(cudaMalloc(&radmin, size));
    check_cudart(cudaMalloc(&epsilon, size));
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
    copyin_data(radmin, radvec.data(), jcount * jcount);
    copyin_data(epsilon, epsvec.data(), jcount * jcount);

    size = n * rs;
    check_cudart(cudaMalloc(&vlam, size));
    std::vector<double> vlamvec(n);
    for (int i = 0; i < n; ++i) {
      if (mutant::mut[i]) {
        vlamvec[i] = mutant::vlambda;
      } else {
        vlamvec[i] = 1;
      }
    }
    copyin_data(vlam, vlamvec.data(), n);

    check_cudart(cudaMalloc(&ev, rs));
    check_cudart(cudaMalloc(&nev, sizeof(int)));
    check_cudart(cudaMalloc(&vir_ev, rs * 9));
  }
}
}
TINKER_NAMESPACE_END
