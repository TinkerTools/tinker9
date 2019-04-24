#ifndef TINKER_MOD_DOMEGA_HH_
#define TINKER_MOD_DOMEGA_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace domega {
extern double*& tesum;
extern double*& teb;
extern double*& tea;
extern double*& teba;
extern double*& teub;
extern double*& teaa;
extern double*& teopb;
extern double*& teopd;
extern double*& teid;
extern double*& teit;
extern double*& tet;
extern double*& tept;
extern double*& tebt;
extern double*& teat;
extern double*& tett;
extern double*& tev;
extern double*& ter;
extern double*& tedsp;
extern double*& tec;
extern double*& tecd;
extern double*& ted;
extern double*& tem;
extern double*& tep;
extern double*& tect;
extern double*& terxf;
extern double*& tes;
extern double*& telf;
extern double*& teg;
extern double*& tex;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(domega, tesum);
extern "C" double* m_tinker_mod(domega, teb);
extern "C" double* m_tinker_mod(domega, tea);
extern "C" double* m_tinker_mod(domega, teba);
extern "C" double* m_tinker_mod(domega, teub);
extern "C" double* m_tinker_mod(domega, teaa);
extern "C" double* m_tinker_mod(domega, teopb);
extern "C" double* m_tinker_mod(domega, teopd);
extern "C" double* m_tinker_mod(domega, teid);
extern "C" double* m_tinker_mod(domega, teit);
extern "C" double* m_tinker_mod(domega, tet);
extern "C" double* m_tinker_mod(domega, tept);
extern "C" double* m_tinker_mod(domega, tebt);
extern "C" double* m_tinker_mod(domega, teat);
extern "C" double* m_tinker_mod(domega, tett);
extern "C" double* m_tinker_mod(domega, tev);
extern "C" double* m_tinker_mod(domega, ter);
extern "C" double* m_tinker_mod(domega, tedsp);
extern "C" double* m_tinker_mod(domega, tec);
extern "C" double* m_tinker_mod(domega, tecd);
extern "C" double* m_tinker_mod(domega, ted);
extern "C" double* m_tinker_mod(domega, tem);
extern "C" double* m_tinker_mod(domega, tep);
extern "C" double* m_tinker_mod(domega, tect);
extern "C" double* m_tinker_mod(domega, terxf);
extern "C" double* m_tinker_mod(domega, tes);
extern "C" double* m_tinker_mod(domega, telf);
extern "C" double* m_tinker_mod(domega, teg);
extern "C" double* m_tinker_mod(domega, tex);

double*& tesum = m_tinker_mod(domega, tesum);
double*& teb = m_tinker_mod(domega, teb);
double*& tea = m_tinker_mod(domega, tea);
double*& teba = m_tinker_mod(domega, teba);
double*& teub = m_tinker_mod(domega, teub);
double*& teaa = m_tinker_mod(domega, teaa);
double*& teopb = m_tinker_mod(domega, teopb);
double*& teopd = m_tinker_mod(domega, teopd);
double*& teid = m_tinker_mod(domega, teid);
double*& teit = m_tinker_mod(domega, teit);
double*& tet = m_tinker_mod(domega, tet);
double*& tept = m_tinker_mod(domega, tept);
double*& tebt = m_tinker_mod(domega, tebt);
double*& teat = m_tinker_mod(domega, teat);
double*& tett = m_tinker_mod(domega, tett);
double*& tev = m_tinker_mod(domega, tev);
double*& ter = m_tinker_mod(domega, ter);
double*& tedsp = m_tinker_mod(domega, tedsp);
double*& tec = m_tinker_mod(domega, tec);
double*& tecd = m_tinker_mod(domega, tecd);
double*& ted = m_tinker_mod(domega, ted);
double*& tem = m_tinker_mod(domega, tem);
double*& tep = m_tinker_mod(domega, tep);
double*& tect = m_tinker_mod(domega, tect);
double*& terxf = m_tinker_mod(domega, terxf);
double*& tes = m_tinker_mod(domega, tes);
double*& telf = m_tinker_mod(domega, telf);
double*& teg = m_tinker_mod(domega, teg);
double*& tex = m_tinker_mod(domega, tex);
#endif

} TINKER_NAMESPACE_END

#endif
