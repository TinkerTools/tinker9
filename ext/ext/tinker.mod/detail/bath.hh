#ifndef TINKER_MOD_BATH_HH_
#define TINKER_MOD_BATH_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace bath {
const int maxnose = 4;
extern int& voltrial;
extern double& kelvin;
extern double& atmsph;
extern double& tautemp;
extern double& taupres;
extern double& compress;
extern double& collide;
extern double& eta;
extern double& volmove;
extern double& vbar;
extern double& qbar;
extern double& gbar;
extern double (&vnh)[maxnose];
extern double (&qnh)[maxnose];
extern double (&gnh)[maxnose];
extern int& isothermal;
extern int& isobaric;
extern int& anisotrop;
extern char (&volscale)[9];
extern char (&barostat)[11];
extern char (&thermostat)[11];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(bath, voltrial);
extern "C" double m_tinker_mod(bath, kelvin);
extern "C" double m_tinker_mod(bath, atmsph);
extern "C" double m_tinker_mod(bath, tautemp);
extern "C" double m_tinker_mod(bath, taupres);
extern "C" double m_tinker_mod(bath, compress);
extern "C" double m_tinker_mod(bath, collide);
extern "C" double m_tinker_mod(bath, eta);
extern "C" double m_tinker_mod(bath, volmove);
extern "C" double m_tinker_mod(bath, vbar);
extern "C" double m_tinker_mod(bath, qbar);
extern "C" double m_tinker_mod(bath, gbar);
extern "C" double m_tinker_mod(bath, vnh)[maxnose];
extern "C" double m_tinker_mod(bath, qnh)[maxnose];
extern "C" double m_tinker_mod(bath, gnh)[maxnose];
extern "C" int m_tinker_mod(bath, isothermal);
extern "C" int m_tinker_mod(bath, isobaric);
extern "C" int m_tinker_mod(bath, anisotrop);
extern "C" char m_tinker_mod(bath, volscale)[9];
extern "C" char m_tinker_mod(bath, barostat)[11];
extern "C" char m_tinker_mod(bath, thermostat)[11];

int& voltrial = m_tinker_mod(bath, voltrial);
double& kelvin = m_tinker_mod(bath, kelvin);
double& atmsph = m_tinker_mod(bath, atmsph);
double& tautemp = m_tinker_mod(bath, tautemp);
double& taupres = m_tinker_mod(bath, taupres);
double& compress = m_tinker_mod(bath, compress);
double& collide = m_tinker_mod(bath, collide);
double& eta = m_tinker_mod(bath, eta);
double& volmove = m_tinker_mod(bath, volmove);
double& vbar = m_tinker_mod(bath, vbar);
double& qbar = m_tinker_mod(bath, qbar);
double& gbar = m_tinker_mod(bath, gbar);
double (&vnh)[maxnose] = m_tinker_mod(bath, vnh);
double (&qnh)[maxnose] = m_tinker_mod(bath, qnh);
double (&gnh)[maxnose] = m_tinker_mod(bath, gnh);
int& isothermal = m_tinker_mod(bath, isothermal);
int& isobaric = m_tinker_mod(bath, isobaric);
int& anisotrop = m_tinker_mod(bath, anisotrop);
char (&volscale)[9] = m_tinker_mod(bath, volscale);
char (&barostat)[11] = m_tinker_mod(bath, barostat);
char (&thermostat)[11] = m_tinker_mod(bath, thermostat);
#endif

} TINKER_NAMESPACE_END

#endif
