#ifndef TINKER_MOD_KATOMS_HH_
#define TINKER_MOD_KATOMS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace katoms {
extern int*& atmcls;
extern int*& atmnum;
extern int*& ligand;
extern double*& weight;
extern char (*&symbol)[3];
extern char (*&describe)[24];

#ifdef TINKER_MOD_CPP_
extern "C" int* m_tinker_mod(katoms, atmcls);
extern "C" int* m_tinker_mod(katoms, atmnum);
extern "C" int* m_tinker_mod(katoms, ligand);
extern "C" double* m_tinker_mod(katoms, weight);
extern "C" char (*m_tinker_mod(katoms, symbol))[3];
extern "C" char (*m_tinker_mod(katoms, describe))[24];

int*& atmcls = m_tinker_mod(katoms, atmcls);
int*& atmnum = m_tinker_mod(katoms, atmnum);
int*& ligand = m_tinker_mod(katoms, ligand);
double*& weight = m_tinker_mod(katoms, weight);
char (*&symbol)[3] = m_tinker_mod(katoms, symbol);
char (*&describe)[24] = m_tinker_mod(katoms, describe);
#endif

} TINKER_NAMESPACE_END

#endif
