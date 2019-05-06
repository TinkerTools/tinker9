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
}
TINKER_NAMESPACE_END

#endif
