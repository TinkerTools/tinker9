#include "acc_seq.h"
#include "util_mdstate.h"

#ifdef TINKER_HOSTONLY
#  include <algorithm>
TINKER_NAMESPACE_BEGIN
static void sort_v1_(int* arr, int len) { std::sort(arr, arr + len); }
TINKER_NAMESPACE_END
#else
TINKER_NAMESPACE_BEGIN
#  define m_swap_(a, b)                                                        \
    {                                                                          \
      auto tmp = b;                                                            \
      b = a;                                                                   \
      a = tmp;                                                                 \
    }
#  define m_max_heap_(arr, start, end)                                         \
    {                                                                          \
      int dad = start;                                                         \
      int son = dad * 2 + 1;                                                   \
      while (son <= end) {                                                     \
        if (son + 1 <= end && arr[son] < arr[son + 1])                         \
          ++son;                                                               \
        if (arr[dad] > arr[son])                                               \
          return;                                                              \
        else {                                                                 \
          m_swap_(arr[dad], arr[son]);                                         \
          dad = son;                                                           \
          son = dad * 2 + 1;                                                   \
        }                                                                      \
      }                                                                        \
    }

#pragma acc routine seq
static void sort_v1_(int* arr, int len) {
  // heapsort
  for (int i = len / 2 - 1; i >= 0; --i) {
    m_max_heap_(arr, i, len - 1);
  }

  for (int i = len - 1; i > 0; --i) {
    m_swap_(arr[0], arr[i]);
    m_max_heap_(arr, 0, i - 1);
  }
}
#  undef m_swap_
#  undef m_max_heap_
TINKER_NAMESPACE_END
#endif

TINKER_NAMESPACE_BEGIN

//======================================================================
// double loop

static void build_double_loop_(nblist_t* lst) {
  #pragma acc parallel loop independent deviceptr(lst)
  for (int i = 0; i < n; ++i) {
    lst->nlst[i] = n - i - 1;
    lst->lst[i] = i + 1;
  }
}

// static void update_double_loop_() {}

//======================================================================
// version 1
// see also nblist.f

static void build_v1_(const nblist_t& st, nblist_t* lst) {
  const int maxnlst = st.maxnlst;
  const real buf2 = REAL_SQ(st.cutoff + st.buffer);

  #pragma acc parallel loop independent deviceptr(lst,box)
  for (int i = 0; i < n; ++i) {
    real xi = lst->x[i];
    real yi = lst->y[i];
    real zi = lst->z[i];
    lst->xold[i] = xi;
    lst->yold[i] = yi;
    lst->zold[i] = zi;

    lst->nlst[i] = 0;
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

static void displace_v1_(const nblist_t& st, nblist_t* lst) {
  const real lbuf2 = REAL_SQ(0.5f * st.buffer);
  #pragma acc parallel loop independent deviceptr(lst,box)
  for (int i = 0; i < n; ++i) {
    real xi = lst->x[i];
    real yi = lst->y[i];
    real zi = lst->z[i];
    real xr = xi - lst->xold[i];
    real yr = yi - lst->yold[i];
    real zr = zi - lst->zold[i];
    imagen(xr, yr, zr, box);
    lst->update[i] = 0;
    real r2 = xr * xr + yr * yr + zr * zr;
    if (r2 >= lbuf2) {
      lst->update[i] = 1;
      lst->xold[i] = xi;
      lst->yold[i] = yi;
      lst->zold[i] = zi;
    }
  }
}

static void update_v1_(const nblist_t& st, nblist_t* lst) {

  // test sites for displacement exceeding half the buffer

  displace_v1_(st, lst);

  // rebuild the higher numbered neighbors for updated sites

  const int maxnlst = st.maxnlst;
  const real buf2 = REAL_SQ(st.cutoff + st.buffer);
  const real bufx = REAL_SQ(st.cutoff + 2 * st.buffer);

  #pragma acc kernels deviceptr(lst,box)
  {
    #pragma acc loop independent
    for (int i = 0; i < n; ++i) {
      if (lst->update[i]) {
        real xi = lst->xold[i];
        real yi = lst->yold[i];
        real zi = lst->zold[i];
        lst->nlst[i] = 0;
        #pragma acc loop seq
        for (int k = i + 1; k < n; ++k) {
          real xr = xi - lst->xold[k];
          real yr = yi - lst->yold[k];
          real zr = zi - lst->zold[k];
          imagen(xr, yr, zr, box);
          real r2 = xr * xr + yr * yr + zr * zr;
          if (r2 <= buf2) {
            int current = lst->nlst[i];
            lst->lst[i * maxnlst + current] = k;
            lst->nlst[i] += 1;
          }
        }
      }
    }

    #pragma acc loop seq
    for (int i = 0; i < n; ++i) {
      if (lst->update[i]) {
        real xi = lst->xold[i];
        real yi = lst->yold[i];
        real zi = lst->zold[i];
        #pragma acc loop seq
        for (int k = 0; k < i - 1; ++k) {
          if (!lst->update[k]) {
            real xr = xi - lst->xold[k];
            real yr = yi - lst->yold[k];
            real zr = zi - lst->zold[k];
            imagen(xr, yr, zr, box);
            real r2 = xr * xr + yr * yr + zr * zr;
            if (r2 <= buf2) {
              for (int j = 0; j < lst->nlst[k]; ++j) {
                if (lst->lst[k * maxnlst + j] == i)
                  goto label_20;
              }
              lst->lst[k * maxnlst + lst->nlst[k]] = i;
              lst->nlst[k] += 1;
            label_20:;
            } else if (r2 <= bufx) {
              for (int j = 0; j < lst->nlst[k]; ++j) {
                if (lst->lst[k * maxnlst + j] == i) {
                  lst->nlst[k] -= 1;
                  lst->lst[k * maxnlst + j] =
                      lst->lst[k * maxnlst + lst->nlst[k]];
                  goto label_30;
                }
              }
            label_30:;
            }
          }
        }
        sort_v1_(&lst->lst[i * maxnlst], lst->nlst[i]);
      }
    }
  }
}

//======================================================================

void nblist_build_acc_impl_(const nblist_t& st, nblist_t* lst) {
  if (st.maxnlst == 1) {
    build_double_loop_(lst);
  } else {
    build_v1_(st, lst);
  }
}

// #define TINKER_DEFAULT_NBLIST_UPDATE_(st, lst) build_v1_(st, lst)
#define TINKER_DEFAULT_NBLIST_UPDATE_(st, lst) update_v1_(st, lst)
void nblist_update_acc_impl_(const nblist_t& st, nblist_t* lst) {
  if (st.maxnlst == 1) {
    // update_double_loop_();
  } else {
    TINKER_DEFAULT_NBLIST_UPDATE_(st, lst);
  }
}
#undef TINKER_DEFAULT_NBLIST_UPDATE_
TINKER_NAMESPACE_END
