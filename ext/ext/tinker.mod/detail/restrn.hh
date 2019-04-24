#ifndef TINKER_MOD_RESTRN_HH_
#define TINKER_MOD_RESTRN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace restrn {
extern int& npfix;
extern int& ndfix;
extern int& nafix;
extern int& ntfix;
extern int& ngfix;
extern int& nchir;
extern int*& ipfix;
extern int*& kpfix;
extern int*& idfix;
extern int*& iafix;
extern int*& itfix;
extern int*& igfix;
extern int*& ichir;
extern double& depth;
extern double& width;
extern double& rwall;
extern double*& xpfix;
extern double*& ypfix;
extern double*& zpfix;
extern double*& pfix;
extern double*& dfix;
extern double*& afix;
extern double*& tfix;
extern double*& gfix;
extern double*& chir;
extern int& use_basin;
extern int& use_wall;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(restrn, npfix);
extern "C" int m_tinker_mod(restrn, ndfix);
extern "C" int m_tinker_mod(restrn, nafix);
extern "C" int m_tinker_mod(restrn, ntfix);
extern "C" int m_tinker_mod(restrn, ngfix);
extern "C" int m_tinker_mod(restrn, nchir);
extern "C" int* m_tinker_mod(restrn, ipfix);
extern "C" int* m_tinker_mod(restrn, kpfix);
extern "C" int* m_tinker_mod(restrn, idfix);
extern "C" int* m_tinker_mod(restrn, iafix);
extern "C" int* m_tinker_mod(restrn, itfix);
extern "C" int* m_tinker_mod(restrn, igfix);
extern "C" int* m_tinker_mod(restrn, ichir);
extern "C" double m_tinker_mod(restrn, depth);
extern "C" double m_tinker_mod(restrn, width);
extern "C" double m_tinker_mod(restrn, rwall);
extern "C" double* m_tinker_mod(restrn, xpfix);
extern "C" double* m_tinker_mod(restrn, ypfix);
extern "C" double* m_tinker_mod(restrn, zpfix);
extern "C" double* m_tinker_mod(restrn, pfix);
extern "C" double* m_tinker_mod(restrn, dfix);
extern "C" double* m_tinker_mod(restrn, afix);
extern "C" double* m_tinker_mod(restrn, tfix);
extern "C" double* m_tinker_mod(restrn, gfix);
extern "C" double* m_tinker_mod(restrn, chir);
extern "C" int m_tinker_mod(restrn, use_basin);
extern "C" int m_tinker_mod(restrn, use_wall);

int& npfix = m_tinker_mod(restrn, npfix);
int& ndfix = m_tinker_mod(restrn, ndfix);
int& nafix = m_tinker_mod(restrn, nafix);
int& ntfix = m_tinker_mod(restrn, ntfix);
int& ngfix = m_tinker_mod(restrn, ngfix);
int& nchir = m_tinker_mod(restrn, nchir);
int*& ipfix = m_tinker_mod(restrn, ipfix);
int*& kpfix = m_tinker_mod(restrn, kpfix);
int*& idfix = m_tinker_mod(restrn, idfix);
int*& iafix = m_tinker_mod(restrn, iafix);
int*& itfix = m_tinker_mod(restrn, itfix);
int*& igfix = m_tinker_mod(restrn, igfix);
int*& ichir = m_tinker_mod(restrn, ichir);
double& depth = m_tinker_mod(restrn, depth);
double& width = m_tinker_mod(restrn, width);
double& rwall = m_tinker_mod(restrn, rwall);
double*& xpfix = m_tinker_mod(restrn, xpfix);
double*& ypfix = m_tinker_mod(restrn, ypfix);
double*& zpfix = m_tinker_mod(restrn, zpfix);
double*& pfix = m_tinker_mod(restrn, pfix);
double*& dfix = m_tinker_mod(restrn, dfix);
double*& afix = m_tinker_mod(restrn, afix);
double*& tfix = m_tinker_mod(restrn, tfix);
double*& gfix = m_tinker_mod(restrn, gfix);
double*& chir = m_tinker_mod(restrn, chir);
int& use_basin = m_tinker_mod(restrn, use_basin);
int& use_wall = m_tinker_mod(restrn, use_wall);
#endif

} TINKER_NAMESPACE_END

#endif
