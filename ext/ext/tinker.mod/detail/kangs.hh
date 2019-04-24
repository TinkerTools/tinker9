#ifndef TINKER_MOD_KANGS_HH_
#define TINKER_MOD_KANGS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kangs {
const int maxna = 2000;
const int maxna5 = 500;
const int maxna4 = 500;
const int maxna3 = 500;
const int maxnaf = 500;
extern double (&acon)[maxna];
extern double (&acon5)[maxna5];
extern double (&acon4)[maxna4];
extern double (&acon3)[maxna3];
extern double (&aconf)[maxnaf];
extern double (&ang)[maxna][3];
extern double (&ang5)[maxna5][3];
extern double (&ang4)[maxna4][3];
extern double (&ang3)[maxna3][3];
extern double (&angf)[maxnaf][2];
extern char (&ka)[maxna][12];
extern char (&ka5)[maxna5][12];
extern char (&ka4)[maxna4][12];
extern char (&ka3)[maxna3][12];
extern char (&kaf)[maxnaf][12];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kangs, acon)[maxna];
extern "C" double m_tinker_mod(kangs, acon5)[maxna5];
extern "C" double m_tinker_mod(kangs, acon4)[maxna4];
extern "C" double m_tinker_mod(kangs, acon3)[maxna3];
extern "C" double m_tinker_mod(kangs, aconf)[maxnaf];
extern "C" double m_tinker_mod(kangs, ang)[maxna][3];
extern "C" double m_tinker_mod(kangs, ang5)[maxna5][3];
extern "C" double m_tinker_mod(kangs, ang4)[maxna4][3];
extern "C" double m_tinker_mod(kangs, ang3)[maxna3][3];
extern "C" double m_tinker_mod(kangs, angf)[maxnaf][2];
extern "C" char m_tinker_mod(kangs, ka)[maxna][12];
extern "C" char m_tinker_mod(kangs, ka5)[maxna5][12];
extern "C" char m_tinker_mod(kangs, ka4)[maxna4][12];
extern "C" char m_tinker_mod(kangs, ka3)[maxna3][12];
extern "C" char m_tinker_mod(kangs, kaf)[maxnaf][12];

double (&acon)[maxna] = m_tinker_mod(kangs, acon);
double (&acon5)[maxna5] = m_tinker_mod(kangs, acon5);
double (&acon4)[maxna4] = m_tinker_mod(kangs, acon4);
double (&acon3)[maxna3] = m_tinker_mod(kangs, acon3);
double (&aconf)[maxnaf] = m_tinker_mod(kangs, aconf);
double (&ang)[maxna][3] = m_tinker_mod(kangs, ang);
double (&ang5)[maxna5][3] = m_tinker_mod(kangs, ang5);
double (&ang4)[maxna4][3] = m_tinker_mod(kangs, ang4);
double (&ang3)[maxna3][3] = m_tinker_mod(kangs, ang3);
double (&angf)[maxnaf][2] = m_tinker_mod(kangs, angf);
char (&ka)[maxna][12] = m_tinker_mod(kangs, ka);
char (&ka5)[maxna5][12] = m_tinker_mod(kangs, ka5);
char (&ka4)[maxna4][12] = m_tinker_mod(kangs, ka4);
char (&ka3)[maxna3][12] = m_tinker_mod(kangs, ka3);
char (&kaf)[maxnaf][12] = m_tinker_mod(kangs, kaf);
#endif

} TINKER_NAMESPACE_END

#endif
