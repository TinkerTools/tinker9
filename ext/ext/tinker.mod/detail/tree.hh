#ifndef TINKER_MOD_TREE_HH_
#define TINKER_MOD_TREE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace tree {
const int maxpss = 500;
extern int& nlevel;
extern double& etree;
extern double (&ilevel)[maxpss+1];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(tree, nlevel);
extern "C" double m_tinker_mod(tree, etree);
extern "C" double m_tinker_mod(tree, ilevel)[maxpss+1];

int& nlevel = m_tinker_mod(tree, nlevel);
double& etree = m_tinker_mod(tree, etree);
double (&ilevel)[maxpss+1] = m_tinker_mod(tree, ilevel);
#endif

} TINKER_NAMESPACE_END

#endif
