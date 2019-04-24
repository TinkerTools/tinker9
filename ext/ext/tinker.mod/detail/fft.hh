#ifndef TINKER_MOD_FFT_HH_
#define TINKER_MOD_FFT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace fft {
const int maxprime = 15;
extern int (&iprime)[3][maxprime];
extern unsigned long long& planf;
extern unsigned long long& planb;
extern double*& ffttable;
extern char (&ffttyp)[7];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(fft, iprime)[3][maxprime];
extern "C" unsigned long long m_tinker_mod(fft, planf);
extern "C" unsigned long long m_tinker_mod(fft, planb);
extern "C" double* m_tinker_mod(fft, ffttable);
extern "C" char m_tinker_mod(fft, ffttyp)[7];

int (&iprime)[3][maxprime] = m_tinker_mod(fft, iprime);
unsigned long long& planf = m_tinker_mod(fft, planf);
unsigned long long& planb = m_tinker_mod(fft, planb);
double*& ffttable = m_tinker_mod(fft, ffttable);
char (&ffttyp)[7] = m_tinker_mod(fft, ffttyp);
#endif

} TINKER_NAMESPACE_END

#endif
