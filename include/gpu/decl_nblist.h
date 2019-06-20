#ifndef TINKER_GPU_DECL_NBLIST_H_
#define TINKER_GPU_DECL_NBLIST_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int list_null = 0;
const int list_double_loop = 1;
const int list_nblist = 2;

struct nblist_st {
  int* nlst;
  int* lst;
  real *xold, *yold, *zold;
  const real* x;
  const real* y;
  const real* z;
  int maxnlst;
  real cutoff, buffer;
};

void nblist_build(const nblist_st&, nblist_st*);
void nblist_update(const nblist_st&, nblist_st*);

int use_vdw_list();
int use_disp_list();
int use_charge_list();
int use_mpole_list();
int use_usolv_list();

extern nblist_st vlist_obj_;
extern nblist_st* vlst;
extern nblist_st dlist_obj_;
extern nblist_st* dlst;
extern nblist_st clist_obj_;
extern nblist_st* clst;
extern nblist_st mlist_obj_;
extern nblist_st* mlst;
extern nblist_st ulist_obj_;
extern nblist_st* ulst;

void nblist_data(int op);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_vlist_update();
void tinker_gpu_mlist_update();
}

#endif
