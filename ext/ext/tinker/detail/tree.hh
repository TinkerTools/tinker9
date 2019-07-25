#ifndef TINKER_MOD_TREE_HH_
#define TINKER_MOD_TREE_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace tree {
const int maxpss = 500;
extern int& nlevel;
extern double& etree;
extern double (&ilevel)[maxpss+1];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(tree, nlevel);
extern "C" double TINKER_MOD(tree, etree);
extern "C" double TINKER_MOD(tree, ilevel)[maxpss+1];

int& nlevel = TINKER_MOD(tree, nlevel);
double& etree = TINKER_MOD(tree, etree);
double (&ilevel)[maxpss+1] = TINKER_MOD(tree, ilevel);
#endif
} TINKER_NAMESPACE_END

#endif
