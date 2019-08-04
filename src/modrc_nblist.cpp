#include "gpu/e_vdw.h"
#include "mod_md.h"
#include "mod_nblist.h"
#include "util_array.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
static int use_vdw_list() {
  int ret = 0;
  if (use_potent(vdw_term))
    ++ret;
  else
    return ret;

  if (limits::use_vlist)
    ++ret;
  return ret;
}

static int use_disp_list() {
  int ret = 0;
  if (potent::use_disp)
    ++ret;
  else
    return ret;

  if (limits::use_dlist)
    ++ret;
  return ret;
}

static int use_charge_list() {
  int ret = 0;
  if (potent::use_charge || potent::use_solv)
    ++ret;
  else
    return ret;

  if (limits::use_clist)
    ++ret;
  return ret;
}

static int use_mpole_list() {
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

static int use_usolv_list() {
  int ret = 0;
  if (potent::use_polar)
    ++ret;
  else
    return ret;

  if (limits::use_ulist)
    ++ret;
  return ret;
}

// see also cutoffs.f
// In the gas phase calculation where neighbor list is not used, we should
// always first check the value of maxn.
// If maxn is equal to 1, it means the value of cutoff can even be INF.
static int nblist_maxlst_(int maxn, double cutoff, double buffer) {
  if (maxn > 1) {
    double buf = (cutoff + buffer);
    int limit = buf * buf * buf + 100;
    int ans = std::min(limit, maxn);
    if (ans > 1) {
      const int magic = 32;
      ans = (ans + magic - 1) / magic;
      ans *= magic;
    }
    return ans;
  } else {
    return 1;
  }
}

static void nblist_op_dealloc_(NBList& st, NBList*& list) {
  dealloc_array(st.nlst);
  dealloc_array(st.lst);
  dealloc_array(st.update);
  dealloc_array(st.xold);
  dealloc_array(st.yold);
  dealloc_array(st.zold);
  dealloc_array(list);
}

static void nblist_op_alloc_(NBList& st, NBList*& list, int maxn, double cutoff,
                             double buffer, const real* _x, const real* _y,
                             const real* _z) {
  const size_t rs = sizeof(int);
  size_t size;

  size = n * rs;
  alloc_array(&st.nlst, size);

  int maxlst = nblist_maxlst_(maxn, cutoff, buffer);
  size = maxlst * n * rs;
  alloc_array(&st.lst, size);

  if (maxlst == 1) {
    st.update = nullptr;
    st.xold = nullptr;
    st.yold = nullptr;
    st.zold = nullptr;
  } else {
    size = n * rs;
    alloc_array(&st.update, size);
    size = n * sizeof(real);
    alloc_array(&st.xold, size);
    alloc_array(&st.yold, size);
    alloc_array(&st.zold, size);
  }

  st.x = _x;
  st.y = _y;
  st.z = _z;

  st.maxnlst = maxlst;
  st.cutoff = cutoff;
  st.buffer = buffer;

  size = sizeof(NBList);
  alloc_array(&list, size);
  copyin_bytes(list, &st, size);
}

extern void nblist_build_acc_impl_(const NBList& st, NBList* lst);
extern void nblist_update_acc_impl_(const NBList&, NBList*);
void nblist_data(rc_op op) {
  int maxnlst = 0;
  int u = 0;

  // vlist
  u = use_vdw_list();
  if (u) {
    if (op & rc_dealloc)
      nblist_op_dealloc_(vlist_obj_, vlst);

    if (op & rc_alloc) {
      maxnlst = 2500;
      if (u == NBList::double_loop)
        maxnlst = 1;
      nblist_op_alloc_(vlist_obj_, vlst, maxnlst, limits::vdwcut,
                       neigh::lbuffer, xred, yred, zred);
    }

    if (op & rc_init) {
      evdw_reduce_xyz();
      nblist_build_acc_impl_(vlist_obj_, vlst);
    }

    if (op & rc_man::evolve) {
      // assuming evdw_reduce_xyz() has been called in the energy routine
      nblist_update_acc_impl_(vlist_obj_, vlst);
    }
  }

  // dlist
  u = use_disp_list();
  if (u) {
    if (op & rc_dealloc)
      nblist_op_dealloc_(dlist_obj_, dlst);

    if (op & rc_alloc) {
      maxnlst = 2500;
      if (u == NBList::double_loop)
        maxnlst = 1;
      nblist_op_alloc_(dlist_obj_, dlst, maxnlst, limits::dispcut,
                       neigh::lbuffer, x, y, z);
    }

    if (op & rc_init) {
    }
  }

  // clist
  u = use_charge_list();
  if (u) {
    if (op & rc_dealloc)
      nblist_op_dealloc_(clist_obj_, clst);

    if (op & rc_alloc) {
      maxnlst = 2500;
      if (u == NBList::double_loop)
        maxnlst = 1;
      nblist_op_alloc_(clist_obj_, clst, maxnlst, limits::chgcut,
                       neigh::lbuffer, x, y, z);
    }

    if (op & rc_init) {
    }
  }

  // mlist
  u = use_mpole_list();
  if (u) {
    if (op & rc_dealloc)
      nblist_op_dealloc_(mlist_obj_, mlst);

    if (op & rc_alloc) {
      maxnlst = 2500;
      if (u == NBList::double_loop)
        maxnlst = 1;
      nblist_op_alloc_(mlist_obj_, mlst, maxnlst, limits::mpolecut,
                       neigh::lbuffer, x, y, z);
    }

    if (op & rc_init)
      nblist_build_acc_impl_(mlist_obj_, mlst);

    if (op & rc_man::evolve) {
      if (use_data & calc::traj) {
        mlist_obj_.x = x;
        mlist_obj_.y = y;
        mlist_obj_.z = z;
        copyin_bytes(mlst, &mlist_obj_, sizeof(NBList));
      }
      nblist_update_acc_impl_(mlist_obj_, mlst);
    }
  }

  // ulist
  u = use_usolv_list();
  if (u) {
    if (op & rc_dealloc)
      nblist_op_dealloc_(ulist_obj_, ulst);

    if (op & rc_alloc) {
      maxnlst = 500;
      if (u == NBList::double_loop)
        maxnlst = 1;
      nblist_op_alloc_(ulist_obj_, ulst, maxnlst, limits::usolvcut,
                       neigh::pbuffer, x, y, z);
    }

    if (op & rc_init) {
    }
  }
}
TINKER_NAMESPACE_END
