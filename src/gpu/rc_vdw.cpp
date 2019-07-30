#include "gpu/e_vdw.h"
#include "util_io.h"
#include "util_mdstate.h"
#include "util_potent.h"
#include "util_switch.h"
#include <ext/tinker/tinker_mod.h>
#include <map>

TINKER_NAMESPACE_BEGIN
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

real ghal, dhal;
real scexp, scalpha;
int vcouple;
real v2scale, v3scale, v4scale, v5scale;

int* ired;
real* kred;
real *xred, *yred, *zred;

int *jvdw, *njvdw;
real *radmin, *epsilon;

real* vlam;

real* ev;
int* nev;
real* vir_ev;

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
  if (!use_potent(vdw_term))
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

    free_nev(nev, ev, vir_ev);
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

    alloc_nev(&nev, &ev, &vir_ev);
  }

  if (rc & rc_copyin) {
    get_evdw_type(vdwtyp);

    switch_cut_off(switch_vdw, vdw_switch_cut, vdw_switch_off);

    ghal = vdwpot::ghal;
    dhal = vdwpot::dhal;
    scexp = mutant::scexp;
    scalpha = mutant::scalpha;
    vcouple = mutant::vcouple;

    v2scale = vdwpot::v2scale;
    v3scale = vdwpot::v3scale;
    v4scale = vdwpot::v4scale;
    v5scale = vdwpot::v5scale;

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

extern void evdw_lj_acc_impl_(int vers);
void evdw_lj(int vers) { evdw_lj_acc_impl_(vers); }

extern void evdw_buck_acc_impl_(int vers);
void evdw_buck(int vers) { evdw_buck_acc_impl_(vers); }

extern void evdw_mm3hb_acc_impl_(int vers);
void evdw_mm3hb(int vers) { evdw_mm3hb_acc_impl_(vers); }

extern void evdw_hal_acc_impl_(int vers);
void evdw_hal(int vers) { evdw_hal_acc_impl_(vers); }

extern void evdw_gauss_acc_impl_(int vers);
void evdw_gauss(int vers) { evdw_gauss_acc_impl_(vers); }

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
TINKER_NAMESPACE_END
