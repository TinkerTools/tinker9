#ifndef TINKER_MOD_NBLIST_H_
#define TINKER_MOD_NBLIST_H_

#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
struct nblist_t {
  enum { null = 0, double_loop = 1, nblist = 2 };

  int* nlst;
  int* lst;
  int* update;
  real *xold, *yold, *zold;
  const real* x;
  const real* y;
  const real* z;
  int maxnlst;
  real cutoff, buffer;
};

TINKER_EXTERN nblist_t vlist_obj_;
TINKER_EXTERN nblist_t* vlst;
TINKER_EXTERN nblist_t dlist_obj_;
TINKER_EXTERN nblist_t* dlst;
TINKER_EXTERN nblist_t clist_obj_;
TINKER_EXTERN nblist_t* clst;
TINKER_EXTERN nblist_t mlist_obj_;
TINKER_EXTERN nblist_t* mlst;
TINKER_EXTERN nblist_t ulist_obj_;
TINKER_EXTERN nblist_t* ulst;

void nblist_data(rc_t rc);
TINKER_NAMESPACE_END

#endif
