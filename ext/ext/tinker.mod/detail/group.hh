#ifndef TINKER_MOD_GROUP_HH_
#define TINKER_MOD_GROUP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace group {
extern int& ngrp;
extern int*& kgrp;
extern int*& grplist;
extern int*& igrp;
extern double*& grpmass;
extern double*& wgrp;
extern int& use_group;
extern int& use_intra;
extern int& use_inter;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(group, ngrp);
extern "C" int* m_tinker_mod(group, kgrp);
extern "C" int* m_tinker_mod(group, grplist);
extern "C" int* m_tinker_mod(group, igrp);
extern "C" double* m_tinker_mod(group, grpmass);
extern "C" double* m_tinker_mod(group, wgrp);
extern "C" int m_tinker_mod(group, use_group);
extern "C" int m_tinker_mod(group, use_intra);
extern "C" int m_tinker_mod(group, use_inter);

int& ngrp = m_tinker_mod(group, ngrp);
int*& kgrp = m_tinker_mod(group, kgrp);
int*& grplist = m_tinker_mod(group, grplist);
int*& igrp = m_tinker_mod(group, igrp);
double*& grpmass = m_tinker_mod(group, grpmass);
double*& wgrp = m_tinker_mod(group, wgrp);
int& use_group = m_tinker_mod(group, use_group);
int& use_intra = m_tinker_mod(group, use_intra);
int& use_inter = m_tinker_mod(group, use_inter);
#endif

} TINKER_NAMESPACE_END

#endif
