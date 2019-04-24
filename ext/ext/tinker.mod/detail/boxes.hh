#ifndef TINKER_MOD_BOXES_HH_
#define TINKER_MOD_BOXES_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace boxes {
extern double& xbox;
extern double& ybox;
extern double& zbox;
extern double& alpha;
extern double& beta;
extern double& gamma;
extern double& xbox2;
extern double& ybox2;
extern double& zbox2;
extern double& box34;
extern double& volbox;
extern double& beta_sin;
extern double& beta_cos;
extern double& gamma_sin;
extern double& gamma_cos;
extern double& beta_term;
extern double& gamma_term;
extern double (&lvec)[3][3];
extern double (&recip)[3][3];
extern int& orthogonal;
extern int& monoclinic;
extern int& triclinic;
extern int& octahedron;
extern char (&spacegrp)[10];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(boxes, xbox);
extern "C" double m_tinker_mod(boxes, ybox);
extern "C" double m_tinker_mod(boxes, zbox);
extern "C" double m_tinker_mod(boxes, alpha);
extern "C" double m_tinker_mod(boxes, beta);
extern "C" double m_tinker_mod(boxes, gamma);
extern "C" double m_tinker_mod(boxes, xbox2);
extern "C" double m_tinker_mod(boxes, ybox2);
extern "C" double m_tinker_mod(boxes, zbox2);
extern "C" double m_tinker_mod(boxes, box34);
extern "C" double m_tinker_mod(boxes, volbox);
extern "C" double m_tinker_mod(boxes, beta_sin);
extern "C" double m_tinker_mod(boxes, beta_cos);
extern "C" double m_tinker_mod(boxes, gamma_sin);
extern "C" double m_tinker_mod(boxes, gamma_cos);
extern "C" double m_tinker_mod(boxes, beta_term);
extern "C" double m_tinker_mod(boxes, gamma_term);
extern "C" double m_tinker_mod(boxes, lvec)[3][3];
extern "C" double m_tinker_mod(boxes, recip)[3][3];
extern "C" int m_tinker_mod(boxes, orthogonal);
extern "C" int m_tinker_mod(boxes, monoclinic);
extern "C" int m_tinker_mod(boxes, triclinic);
extern "C" int m_tinker_mod(boxes, octahedron);
extern "C" char m_tinker_mod(boxes, spacegrp)[10];

double& xbox = m_tinker_mod(boxes, xbox);
double& ybox = m_tinker_mod(boxes, ybox);
double& zbox = m_tinker_mod(boxes, zbox);
double& alpha = m_tinker_mod(boxes, alpha);
double& beta = m_tinker_mod(boxes, beta);
double& gamma = m_tinker_mod(boxes, gamma);
double& xbox2 = m_tinker_mod(boxes, xbox2);
double& ybox2 = m_tinker_mod(boxes, ybox2);
double& zbox2 = m_tinker_mod(boxes, zbox2);
double& box34 = m_tinker_mod(boxes, box34);
double& volbox = m_tinker_mod(boxes, volbox);
double& beta_sin = m_tinker_mod(boxes, beta_sin);
double& beta_cos = m_tinker_mod(boxes, beta_cos);
double& gamma_sin = m_tinker_mod(boxes, gamma_sin);
double& gamma_cos = m_tinker_mod(boxes, gamma_cos);
double& beta_term = m_tinker_mod(boxes, beta_term);
double& gamma_term = m_tinker_mod(boxes, gamma_term);
double (&lvec)[3][3] = m_tinker_mod(boxes, lvec);
double (&recip)[3][3] = m_tinker_mod(boxes, recip);
int& orthogonal = m_tinker_mod(boxes, orthogonal);
int& monoclinic = m_tinker_mod(boxes, monoclinic);
int& triclinic = m_tinker_mod(boxes, triclinic);
int& octahedron = m_tinker_mod(boxes, octahedron);
char (&spacegrp)[10] = m_tinker_mod(boxes, spacegrp);
#endif

} TINKER_NAMESPACE_END

#endif
