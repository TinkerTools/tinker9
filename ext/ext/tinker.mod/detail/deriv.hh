#ifndef TINKER_MOD_DERIV_HH_
#define TINKER_MOD_DERIV_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace deriv {
extern double*& desum;
extern double*& deb;
extern double*& dea;
extern double*& deba;
extern double*& deub;
extern double*& deaa;
extern double*& deopb;
extern double*& deopd;
extern double*& deid;
extern double*& deit;
extern double*& det;
extern double*& dept;
extern double*& debt;
extern double*& deat;
extern double*& dett;
extern double*& dev;
extern double*& der;
extern double*& dedsp;
extern double*& dec;
extern double*& decd;
extern double*& ded;
extern double*& dem;
extern double*& dep;
extern double*& dect;
extern double*& derxf;
extern double*& des;
extern double*& delf;
extern double*& deg;
extern double*& dex;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(deriv, desum);
extern "C" double* m_tinker_mod(deriv, deb);
extern "C" double* m_tinker_mod(deriv, dea);
extern "C" double* m_tinker_mod(deriv, deba);
extern "C" double* m_tinker_mod(deriv, deub);
extern "C" double* m_tinker_mod(deriv, deaa);
extern "C" double* m_tinker_mod(deriv, deopb);
extern "C" double* m_tinker_mod(deriv, deopd);
extern "C" double* m_tinker_mod(deriv, deid);
extern "C" double* m_tinker_mod(deriv, deit);
extern "C" double* m_tinker_mod(deriv, det);
extern "C" double* m_tinker_mod(deriv, dept);
extern "C" double* m_tinker_mod(deriv, debt);
extern "C" double* m_tinker_mod(deriv, deat);
extern "C" double* m_tinker_mod(deriv, dett);
extern "C" double* m_tinker_mod(deriv, dev);
extern "C" double* m_tinker_mod(deriv, der);
extern "C" double* m_tinker_mod(deriv, dedsp);
extern "C" double* m_tinker_mod(deriv, dec);
extern "C" double* m_tinker_mod(deriv, decd);
extern "C" double* m_tinker_mod(deriv, ded);
extern "C" double* m_tinker_mod(deriv, dem);
extern "C" double* m_tinker_mod(deriv, dep);
extern "C" double* m_tinker_mod(deriv, dect);
extern "C" double* m_tinker_mod(deriv, derxf);
extern "C" double* m_tinker_mod(deriv, des);
extern "C" double* m_tinker_mod(deriv, delf);
extern "C" double* m_tinker_mod(deriv, deg);
extern "C" double* m_tinker_mod(deriv, dex);

double*& desum = m_tinker_mod(deriv, desum);
double*& deb = m_tinker_mod(deriv, deb);
double*& dea = m_tinker_mod(deriv, dea);
double*& deba = m_tinker_mod(deriv, deba);
double*& deub = m_tinker_mod(deriv, deub);
double*& deaa = m_tinker_mod(deriv, deaa);
double*& deopb = m_tinker_mod(deriv, deopb);
double*& deopd = m_tinker_mod(deriv, deopd);
double*& deid = m_tinker_mod(deriv, deid);
double*& deit = m_tinker_mod(deriv, deit);
double*& det = m_tinker_mod(deriv, det);
double*& dept = m_tinker_mod(deriv, dept);
double*& debt = m_tinker_mod(deriv, debt);
double*& deat = m_tinker_mod(deriv, deat);
double*& dett = m_tinker_mod(deriv, dett);
double*& dev = m_tinker_mod(deriv, dev);
double*& der = m_tinker_mod(deriv, der);
double*& dedsp = m_tinker_mod(deriv, dedsp);
double*& dec = m_tinker_mod(deriv, dec);
double*& decd = m_tinker_mod(deriv, decd);
double*& ded = m_tinker_mod(deriv, ded);
double*& dem = m_tinker_mod(deriv, dem);
double*& dep = m_tinker_mod(deriv, dep);
double*& dect = m_tinker_mod(deriv, dect);
double*& derxf = m_tinker_mod(deriv, derxf);
double*& des = m_tinker_mod(deriv, des);
double*& delf = m_tinker_mod(deriv, delf);
double*& deg = m_tinker_mod(deriv, deg);
double*& dex = m_tinker_mod(deriv, dex);
#endif

} TINKER_NAMESPACE_END

#endif
