#pragma once
#include "gen_unit.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
struct KDTree
{
   ~KDTree();
};
using KDTreeUnit = GenericUnit<KDTree, GenericUnitVersion::EnableOnDevice>;
TINKER_EXTERN KDTreeUnit vtree_unit;


void spatial_data(rc_op);
TINKER_NAMESPACE_END
