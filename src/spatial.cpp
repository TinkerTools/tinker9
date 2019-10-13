#include "spatial.h"


TINKER_NAMESPACE_BEGIN
KDTree::~KDTree() {}


extern void spatial_data_cu(rc_op);
void spatial_data(rc_op op)
{
#if TINKER_CUDART
   spatial_data_cu(op);
#endif
}
TINKER_NAMESPACE_END
