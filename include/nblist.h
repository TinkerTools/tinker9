#ifndef TINKER_NBLIST_H_
#define TINKER_NBLIST_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
struct NBList {
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

TINKER_EXTERN NBList vlist_obj_;
TINKER_EXTERN NBList* vlst;
TINKER_EXTERN NBList dlist_obj_;
TINKER_EXTERN NBList* dlst;
TINKER_EXTERN NBList clist_obj_;
TINKER_EXTERN NBList* clst;
TINKER_EXTERN NBList mlist_obj_;
TINKER_EXTERN NBList* mlst;
TINKER_EXTERN NBList ulist_obj_;
TINKER_EXTERN NBList* ulst;

void nblist_data(rc_op op);
TINKER_NAMESPACE_END

#endif
