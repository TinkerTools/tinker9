#ifndef TINKER_MOD_KDIPOL_HH_
#define TINKER_MOD_KDIPOL_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kdipol {
const int maxnd = 1000;
const int maxnd5 = 500;
const int maxnd4 = 500;
const int maxnd3 = 500;
extern double (&dpl)[maxnd];
extern double (&dpl5)[maxnd5];
extern double (&dpl4)[maxnd4];
extern double (&dpl3)[maxnd3];
extern double (&pos)[maxnd];
extern double (&pos5)[maxnd5];
extern double (&pos4)[maxnd4];
extern double (&pos3)[maxnd3];
extern char (&kd)[maxnd][8];
extern char (&kd5)[maxnd5][8];
extern char (&kd4)[maxnd4][8];
extern char (&kd3)[maxnd3][8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kdipol, dpl)[maxnd];
extern "C" double m_tinker_mod(kdipol, dpl5)[maxnd5];
extern "C" double m_tinker_mod(kdipol, dpl4)[maxnd4];
extern "C" double m_tinker_mod(kdipol, dpl3)[maxnd3];
extern "C" double m_tinker_mod(kdipol, pos)[maxnd];
extern "C" double m_tinker_mod(kdipol, pos5)[maxnd5];
extern "C" double m_tinker_mod(kdipol, pos4)[maxnd4];
extern "C" double m_tinker_mod(kdipol, pos3)[maxnd3];
extern "C" char m_tinker_mod(kdipol, kd)[maxnd][8];
extern "C" char m_tinker_mod(kdipol, kd5)[maxnd5][8];
extern "C" char m_tinker_mod(kdipol, kd4)[maxnd4][8];
extern "C" char m_tinker_mod(kdipol, kd3)[maxnd3][8];

double (&dpl)[maxnd] = m_tinker_mod(kdipol, dpl);
double (&dpl5)[maxnd5] = m_tinker_mod(kdipol, dpl5);
double (&dpl4)[maxnd4] = m_tinker_mod(kdipol, dpl4);
double (&dpl3)[maxnd3] = m_tinker_mod(kdipol, dpl3);
double (&pos)[maxnd] = m_tinker_mod(kdipol, pos);
double (&pos5)[maxnd5] = m_tinker_mod(kdipol, pos5);
double (&pos4)[maxnd4] = m_tinker_mod(kdipol, pos4);
double (&pos3)[maxnd3] = m_tinker_mod(kdipol, pos3);
char (&kd)[maxnd][8] = m_tinker_mod(kdipol, kd);
char (&kd5)[maxnd5][8] = m_tinker_mod(kdipol, kd5);
char (&kd4)[maxnd4][8] = m_tinker_mod(kdipol, kd4);
char (&kd3)[maxnd3][8] = m_tinker_mod(kdipol, kd3);
#endif

} TINKER_NAMESPACE_END

#endif
