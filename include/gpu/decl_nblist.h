#ifndef TINKER_GPU_DECL_NBLIST_H_
#define TINKER_GPU_DECL_NBLIST_H_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
typedef struct nblist_def_st_ {
  enum { null = 0, double_loop = 1, nblist = 2 };

  int* nlst;
  int* lst;
  real *xold, *yold, *zold;
  const real* x;
  const real* y;
  const real* z;
  int maxnlst;
  real cutoff, buffer;
} nblist_t;

void nblist_build(const nblist_t&, nblist_t*);
void nblist_update(const nblist_t&, nblist_t*);
void nblist_update_vdw_list();

int use_vdw_list();
int use_disp_list();
int use_charge_list();
int use_mpole_list();
int use_usolv_list();

extern nblist_t vlist_obj_;
extern nblist_t* vlst;
extern nblist_t dlist_obj_;
extern nblist_t* dlst;
extern nblist_t clist_obj_;
extern nblist_t* clst;
extern nblist_t mlist_obj_;
extern nblist_t* mlst;
extern nblist_t ulist_obj_;
extern nblist_t* ulst;

void nblist_data(rc_t rc);
}
TINKER_NAMESPACE_END

#endif
