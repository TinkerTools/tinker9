#ifndef TINKER_MOD_DMA_HH_
#define TINKER_MOD_DMA_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace dma {
extern double*& mp;
extern double*& dpx;
extern double*& dpy;
extern double*& dpz;
extern double*& q20;
extern double*& q21c;
extern double*& q21s;
extern double*& q22c;
extern double*& q22s;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(dma, mp);
extern "C" double* m_tinker_mod(dma, dpx);
extern "C" double* m_tinker_mod(dma, dpy);
extern "C" double* m_tinker_mod(dma, dpz);
extern "C" double* m_tinker_mod(dma, q20);
extern "C" double* m_tinker_mod(dma, q21c);
extern "C" double* m_tinker_mod(dma, q21s);
extern "C" double* m_tinker_mod(dma, q22c);
extern "C" double* m_tinker_mod(dma, q22s);

double*& mp = m_tinker_mod(dma, mp);
double*& dpx = m_tinker_mod(dma, dpx);
double*& dpy = m_tinker_mod(dma, dpy);
double*& dpz = m_tinker_mod(dma, dpz);
double*& q20 = m_tinker_mod(dma, q20);
double*& q21c = m_tinker_mod(dma, q21c);
double*& q21s = m_tinker_mod(dma, q21s);
double*& q22c = m_tinker_mod(dma, q22c);
double*& q22s = m_tinker_mod(dma, q22s);
#endif

} TINKER_NAMESPACE_END

#endif
