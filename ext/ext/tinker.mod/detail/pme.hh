#ifndef TINKER_MOD_PME_HH_
#define TINKER_MOD_PME_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace pme {
extern int& nfft1;
extern int& nfft2;
extern int& nfft3;
extern int& nefft1;
extern int& nefft2;
extern int& nefft3;
extern int& ndfft1;
extern int& ndfft2;
extern int& ndfft3;
extern int& bsorder;
extern int& bseorder;
extern int& bsporder;
extern int& bsdorder;
extern int*& igrid;
extern double*& bsmod1;
extern double*& bsmod2;
extern double*& bsmod3;
extern double*& bsbuild;
extern double*& thetai1;
extern double*& thetai2;
extern double*& thetai3;
extern double*& qgrid;
extern double*& qfac;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(pme, nfft1);
extern "C" int m_tinker_mod(pme, nfft2);
extern "C" int m_tinker_mod(pme, nfft3);
extern "C" int m_tinker_mod(pme, nefft1);
extern "C" int m_tinker_mod(pme, nefft2);
extern "C" int m_tinker_mod(pme, nefft3);
extern "C" int m_tinker_mod(pme, ndfft1);
extern "C" int m_tinker_mod(pme, ndfft2);
extern "C" int m_tinker_mod(pme, ndfft3);
extern "C" int m_tinker_mod(pme, bsorder);
extern "C" int m_tinker_mod(pme, bseorder);
extern "C" int m_tinker_mod(pme, bsporder);
extern "C" int m_tinker_mod(pme, bsdorder);
extern "C" int* m_tinker_mod(pme, igrid);
extern "C" double* m_tinker_mod(pme, bsmod1);
extern "C" double* m_tinker_mod(pme, bsmod2);
extern "C" double* m_tinker_mod(pme, bsmod3);
extern "C" double* m_tinker_mod(pme, bsbuild);
extern "C" double* m_tinker_mod(pme, thetai1);
extern "C" double* m_tinker_mod(pme, thetai2);
extern "C" double* m_tinker_mod(pme, thetai3);
extern "C" double* m_tinker_mod(pme, qgrid);
extern "C" double* m_tinker_mod(pme, qfac);

int& nfft1 = m_tinker_mod(pme, nfft1);
int& nfft2 = m_tinker_mod(pme, nfft2);
int& nfft3 = m_tinker_mod(pme, nfft3);
int& nefft1 = m_tinker_mod(pme, nefft1);
int& nefft2 = m_tinker_mod(pme, nefft2);
int& nefft3 = m_tinker_mod(pme, nefft3);
int& ndfft1 = m_tinker_mod(pme, ndfft1);
int& ndfft2 = m_tinker_mod(pme, ndfft2);
int& ndfft3 = m_tinker_mod(pme, ndfft3);
int& bsorder = m_tinker_mod(pme, bsorder);
int& bseorder = m_tinker_mod(pme, bseorder);
int& bsporder = m_tinker_mod(pme, bsporder);
int& bsdorder = m_tinker_mod(pme, bsdorder);
int*& igrid = m_tinker_mod(pme, igrid);
double*& bsmod1 = m_tinker_mod(pme, bsmod1);
double*& bsmod2 = m_tinker_mod(pme, bsmod2);
double*& bsmod3 = m_tinker_mod(pme, bsmod3);
double*& bsbuild = m_tinker_mod(pme, bsbuild);
double*& thetai1 = m_tinker_mod(pme, thetai1);
double*& thetai2 = m_tinker_mod(pme, thetai2);
double*& thetai3 = m_tinker_mod(pme, thetai3);
double*& qgrid = m_tinker_mod(pme, qgrid);
double*& qfac = m_tinker_mod(pme, qfac);
#endif

} TINKER_NAMESPACE_END

#endif
