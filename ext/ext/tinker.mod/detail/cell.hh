#ifndef TINKER_MOD_CELL_HH_
#define TINKER_MOD_CELL_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace cell {
extern int& ncell;
extern int*& icell;
extern double& xcell;
extern double& ycell;
extern double& zcell;
extern double& xcell2;
extern double& ycell2;
extern double& zcell2;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(cell, ncell);
extern "C" int* m_tinker_mod(cell, icell);
extern "C" double m_tinker_mod(cell, xcell);
extern "C" double m_tinker_mod(cell, ycell);
extern "C" double m_tinker_mod(cell, zcell);
extern "C" double m_tinker_mod(cell, xcell2);
extern "C" double m_tinker_mod(cell, ycell2);
extern "C" double m_tinker_mod(cell, zcell2);

int& ncell = m_tinker_mod(cell, ncell);
int*& icell = m_tinker_mod(cell, icell);
double& xcell = m_tinker_mod(cell, xcell);
double& ycell = m_tinker_mod(cell, ycell);
double& zcell = m_tinker_mod(cell, zcell);
double& xcell2 = m_tinker_mod(cell, xcell2);
double& ycell2 = m_tinker_mod(cell, ycell2);
double& zcell2 = m_tinker_mod(cell, zcell2);
#endif

} TINKER_NAMESPACE_END

#endif
