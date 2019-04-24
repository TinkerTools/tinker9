#ifndef TINKER_MOD_DSPPOT_HH_
#define TINKER_MOD_DSPPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace dsppot {
extern double& dsp2scale;
extern double& dsp3scale;
extern double& dsp4scale;
extern double& dsp5scale;
extern int& use_dcorr;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(dsppot, dsp2scale);
extern "C" double m_tinker_mod(dsppot, dsp3scale);
extern "C" double m_tinker_mod(dsppot, dsp4scale);
extern "C" double m_tinker_mod(dsppot, dsp5scale);
extern "C" int m_tinker_mod(dsppot, use_dcorr);

double& dsp2scale = m_tinker_mod(dsppot, dsp2scale);
double& dsp3scale = m_tinker_mod(dsppot, dsp3scale);
double& dsp4scale = m_tinker_mod(dsppot, dsp4scale);
double& dsp5scale = m_tinker_mod(dsppot, dsp5scale);
int& use_dcorr = m_tinker_mod(dsppot, use_dcorr);
#endif

} TINKER_NAMESPACE_END

#endif
