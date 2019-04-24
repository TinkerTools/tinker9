#ifndef TINKER_MOD_WARP_HH_
#define TINKER_MOD_WARP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace warp {
extern double& deform;
extern double& difft;
extern double& diffv;
extern double& diffc;
extern double*& m2;
extern int& use_smooth;
extern int& use_dem;
extern int& use_gda;
extern int& use_tophat;
extern int& use_stophat;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(warp, deform);
extern "C" double m_tinker_mod(warp, difft);
extern "C" double m_tinker_mod(warp, diffv);
extern "C" double m_tinker_mod(warp, diffc);
extern "C" double* m_tinker_mod(warp, m2);
extern "C" int m_tinker_mod(warp, use_smooth);
extern "C" int m_tinker_mod(warp, use_dem);
extern "C" int m_tinker_mod(warp, use_gda);
extern "C" int m_tinker_mod(warp, use_tophat);
extern "C" int m_tinker_mod(warp, use_stophat);

double& deform = m_tinker_mod(warp, deform);
double& difft = m_tinker_mod(warp, difft);
double& diffv = m_tinker_mod(warp, diffv);
double& diffc = m_tinker_mod(warp, diffc);
double*& m2 = m_tinker_mod(warp, m2);
int& use_smooth = m_tinker_mod(warp, use_smooth);
int& use_dem = m_tinker_mod(warp, use_dem);
int& use_gda = m_tinker_mod(warp, use_gda);
int& use_tophat = m_tinker_mod(warp, use_tophat);
int& use_stophat = m_tinker_mod(warp, use_stophat);
#endif

} TINKER_NAMESPACE_END

#endif
