#ifndef TINKER_GPU_DECL_NBLIST_H_
#define TINKER_GPU_DECL_NBLIST_H_

#include "mod_nblist.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int use_vdw_list();
int use_disp_list();
int use_charge_list();
int use_mpole_list();
int use_usolv_list();

void nblist_data(rc_t rc);
}
TINKER_NAMESPACE_END

#endif
