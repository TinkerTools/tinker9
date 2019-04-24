#ifndef TINKER_MOD_KBONDS_HH_
#define TINKER_MOD_KBONDS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kbonds {
const int maxnb = 2000;
const int maxnb5 = 500;
const int maxnb4 = 500;
const int maxnb3 = 500;
const int maxnel = 500;
extern double (&bcon)[maxnb];
extern double (&bcon5)[maxnb5];
extern double (&bcon4)[maxnb4];
extern double (&bcon3)[maxnb3];
extern double (&blen)[maxnb];
extern double (&blen5)[maxnb5];
extern double (&blen4)[maxnb4];
extern double (&blen3)[maxnb3];
extern double (&dlen)[maxnel];
extern char (&kb)[maxnb][8];
extern char (&kb5)[maxnb5][8];
extern char (&kb4)[maxnb4][8];
extern char (&kb3)[maxnb3][8];
extern char (&kel)[maxnel][12];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kbonds, bcon)[maxnb];
extern "C" double m_tinker_mod(kbonds, bcon5)[maxnb5];
extern "C" double m_tinker_mod(kbonds, bcon4)[maxnb4];
extern "C" double m_tinker_mod(kbonds, bcon3)[maxnb3];
extern "C" double m_tinker_mod(kbonds, blen)[maxnb];
extern "C" double m_tinker_mod(kbonds, blen5)[maxnb5];
extern "C" double m_tinker_mod(kbonds, blen4)[maxnb4];
extern "C" double m_tinker_mod(kbonds, blen3)[maxnb3];
extern "C" double m_tinker_mod(kbonds, dlen)[maxnel];
extern "C" char m_tinker_mod(kbonds, kb)[maxnb][8];
extern "C" char m_tinker_mod(kbonds, kb5)[maxnb5][8];
extern "C" char m_tinker_mod(kbonds, kb4)[maxnb4][8];
extern "C" char m_tinker_mod(kbonds, kb3)[maxnb3][8];
extern "C" char m_tinker_mod(kbonds, kel)[maxnel][12];

double (&bcon)[maxnb] = m_tinker_mod(kbonds, bcon);
double (&bcon5)[maxnb5] = m_tinker_mod(kbonds, bcon5);
double (&bcon4)[maxnb4] = m_tinker_mod(kbonds, bcon4);
double (&bcon3)[maxnb3] = m_tinker_mod(kbonds, bcon3);
double (&blen)[maxnb] = m_tinker_mod(kbonds, blen);
double (&blen5)[maxnb5] = m_tinker_mod(kbonds, blen5);
double (&blen4)[maxnb4] = m_tinker_mod(kbonds, blen4);
double (&blen3)[maxnb3] = m_tinker_mod(kbonds, blen3);
double (&dlen)[maxnel] = m_tinker_mod(kbonds, dlen);
char (&kb)[maxnb][8] = m_tinker_mod(kbonds, kb);
char (&kb5)[maxnb5][8] = m_tinker_mod(kbonds, kb5);
char (&kb4)[maxnb4][8] = m_tinker_mod(kbonds, kb4);
char (&kb3)[maxnb3][8] = m_tinker_mod(kbonds, kb3);
char (&kel)[maxnel][12] = m_tinker_mod(kbonds, kel);
#endif

} TINKER_NAMESPACE_END

#endif
