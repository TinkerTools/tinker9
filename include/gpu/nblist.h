#ifndef TINKER_GPU_NBLIST_H_
#define TINKER_GPU_NBLIST_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int list_none = 0;
const int list_double_loop = 1;
const int list_nblist = 2;

int use_vdw_list();
int use_disp_list();
int use_charge_list();
int use_mpole_list();
int use_usolv_list();

void nblist_construct(const nblist_st&, nblist_st*);
void nblist_update(const nblist_st&, nblist_st*);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_vlist_build();
void tinker_gpu_vlist_update();
}

#endif
