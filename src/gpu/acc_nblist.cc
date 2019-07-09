#include "gpu/acc.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {

//======================================================================
// double loop
static void build_double_loop_(nblist_t* lst) {
  #pragma acc data deviceptr(lst)
  #pragma acc parallel loop
  for (int i = 0; i < n; ++i) {
    lst->nlst[i] = n - i - 1;
    lst->lst[i] = i + 1;
  }
}

// static void update_double_loop_() {}

//======================================================================
// version 1
// see also nblist.f

static void build_v1_(nblist_t* lst) {
  #pragma acc data deviceptr(lst,box)
  #pragma acc parallel loop
  for (int i = 0; i < n; ++i) {
    real xi = lst->x[i];
    real yi = lst->y[i];
    real zi = lst->z[i];
    lst->xold[i] = xi;
    lst->yold[i] = yi;
    lst->zold[i] = zi;

    lst->nlst[i] = 0;
    const int maxnlst = lst->maxnlst;
    const real buf2 = (lst->cutoff + lst->buffer) * (lst->cutoff + lst->buffer);
    #pragma acc loop seq
    for (int k = i + 1; k < n; ++k) {
      real xr = xi - lst->x[k];
      real yr = yi - lst->y[k];
      real zr = zi - lst->z[k];
      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= buf2) {
        const int j = lst->nlst[i];
        lst->nlst[i] += 1;
        lst->lst[i * maxnlst + j] = k;
      }
    }
  }
}

static void update_v1_(nblist_t*) {}

//======================================================================

void nblist_build_acc_impl_(const nblist_t& st, nblist_t* lst) {
  if (st.maxnlst == 1) {
    build_double_loop_(lst);
  } else {
    build_v1_(lst);
  }
}

void nblist_update_acc_impl_(const nblist_t& st, nblist_t* lst) {
  if (st.maxnlst == 1) {
    // update_double_loop_();
  } else {
    update_v1_(lst);
  }
}
}
TINKER_NAMESPACE_END
