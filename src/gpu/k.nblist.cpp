#include "gpu/e.vdw.h"
#include "gpu/mdstate.h"
#include "gpu/nblist.h"
#include "tinker.mod.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int use_vdw_list() {
  int ret = 0;
  if (use_evdw())
    ++ret;
  else
    return ret;

  if (limits::use_vlist)
    ++ret;
  return ret;
}

int use_disp_list() {
  int ret = 0;
  if (potent::use_disp)
    ++ret;
  else
    return ret;

  if (limits::use_dlist)
    ++ret;
  return ret;
}

int use_charge_list() {
  int ret = 0;
  if (potent::use_charge || potent::use_solv)
    ++ret;
  else
    return ret;

  if (limits::use_clist)
    ++ret;
  return ret;
}

int use_mpole_list() {
  int ret = 0;
  if (potent::use_mpole || potent::use_polar || potent::use_chgtrn ||
      potent::use_solv)
    ++ret;
  else
    return ret;

  if (limits::use_mlist)
    ++ret;
  return ret;
}

int use_usolv_list() {
  int ret = 0;
  if (potent::use_polar)
    ++ret;
  else
    return ret;

  if (limits::use_ulist)
    ++ret;
  return ret;
}

void nblist_data_0_(nblist_st& st, nblist_st*& list) {
  check_cudart(cudaFree(st.nlst));
  check_cudart(cudaFree(st.lst));
  check_cudart(cudaFree(st.xold));
  check_cudart(cudaFree(st.yold));
  check_cudart(cudaFree(st.zold));
  check_cudart(cudaFree(list));
}

// see also cutoffs.f
int nblist_maxlst_(int maxn, double cutoff, double buffer) {
  double buf = (cutoff + buffer);
  int limit = buf * buf * buf + 100;
  int ans = std::min(limit, maxn);

  if (ans > 1) {
    const int magic = 32;
    ans = (ans + magic - 1) / magic;
    ans *= magic;
  }
  return ans;
}

void nblist_data_1_(nblist_st& st, nblist_st*& list, int maxn, double cutoff,
                    double buffer, const real* _x, const real* _y,
                    const real* _z) {
  const size_t rs = sizeof(int);
  size_t size;

  size = n * rs;
  check_cudart(cudaMalloc(&st.nlst, size));

  int maxlst = nblist_maxlst_(maxn, cutoff, buffer);
  size = maxlst * n * rs;
  check_cudart(cudaMalloc(&st.lst, size));

  st.maxnlst = maxlst;
  st.cutoff = cutoff;
  st.buffer = buffer;

  if (maxlst == 1) {
    st.xold = nullptr;
    st.yold = nullptr;
    st.zold = nullptr;
  } else {
    size = n * sizeof(real);
    check_cudart(cudaMalloc(&st.xold, size));
    check_cudart(cudaMalloc(&st.yold, size));
    check_cudart(cudaMalloc(&st.zold, size));
  }

  st.x = _x;
  st.y = _y;
  st.z = _z;

  size = sizeof(nblist_st);
  check_cudart(cudaMalloc(&list, size));
  check_cudart(cudaMemcpy(list, &st, size, cudaMemcpyHostToDevice));
}

void nblist_data(int op) {
  int maxnlst = 0;
  int u = 0;

  // vlist
  u = use_vdw_list();
  if (u) {
    if (op == op_destroy)
      nblist_data_0_(vlist_obj_, vlst);

    if (op == op_create) {
      maxnlst = 2500;
      if (u == list_double_loop)
        maxnlst = 1;
      nblist_data_1_(vlist_obj_, vlst, maxnlst, limits::vdwcut, neigh::lbuffer,
                     xred, yred, zred);
    }
  }

  // dlist
  u = use_disp_list();
  if (u) {
    if (op == op_destroy)
      nblist_data_0_(dlist_obj_, dlst);

    if (op == op_create) {
      maxnlst = 2500;
      if (u == list_double_loop)
        maxnlst = 1;
      nblist_data_1_(dlist_obj_, dlst, maxnlst, limits::dispcut, neigh::lbuffer,
                     x, y, z);
    }
  }

  // clist
  u = use_charge_list();
  if (u) {
    if (op == op_destroy)
      nblist_data_0_(clist_obj_, clst);

    if (op == op_create) {
      maxnlst = 2500;
      if (u == list_double_loop)
        maxnlst = 1;
      nblist_data_1_(clist_obj_, clst, maxnlst, limits::chgcut, neigh::lbuffer,
                     x, y, z);
    }
  }

  // mlist
  u = use_mpole_list();
  if (u) {
    if (op == op_destroy)
      nblist_data_0_(mlist_obj_, mlst);

    if (op == op_create) {
      maxnlst = 2500;
      if (u == list_double_loop)
        maxnlst = 1;
      nblist_data_1_(mlist_obj_, mlst, maxnlst, limits::mpolecut,
                     neigh::lbuffer, x, y, z);
    }
  }

  // ulist
  u = use_usolv_list();
  if (u) {
    if (op == op_destroy)
      nblist_data_0_(ulist_obj_, ulst);

    if (op == op_create) {
      maxnlst = 500;
      if (u == list_double_loop)
        maxnlst = 1;
      nblist_data_1_(ulist_obj_, ulst, maxnlst, limits::usolvcut,
                     neigh::pbuffer, x, y, z);
    }
  }
}
}
TINKER_NAMESPACE_END
