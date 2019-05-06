#ifndef TINKER_GPU_MDSTATE_H_
#define TINKER_GPU_MDSTATE_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern int use_data;
extern int n;

extern real* esum;
extern real* vir;
extern real* mass;

extern real *x, *y, *z;
extern real *vx, *vy, *vz;
extern real *ax, *ay, *az;
extern real *gx, *gy, *gz;

extern box_st* box;
extern couple_st couple_obj_;
extern couple_st* couple;

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
}
TINKER_NAMESPACE_END

#endif
