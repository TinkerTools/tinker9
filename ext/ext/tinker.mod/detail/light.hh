#ifndef TINKER_MOD_LIGHT_HH_
#define TINKER_MOD_LIGHT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace light {
extern int& nlight;
extern int*& kbx;
extern int*& kby;
extern int*& kbz;
extern int*& kex;
extern int*& key;
extern int*& kez;
extern int*& locx;
extern int*& locy;
extern int*& locz;
extern int*& rgx;
extern int*& rgy;
extern int*& rgz;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(light, nlight);
extern "C" int* m_tinker_mod(light, kbx);
extern "C" int* m_tinker_mod(light, kby);
extern "C" int* m_tinker_mod(light, kbz);
extern "C" int* m_tinker_mod(light, kex);
extern "C" int* m_tinker_mod(light, key);
extern "C" int* m_tinker_mod(light, kez);
extern "C" int* m_tinker_mod(light, locx);
extern "C" int* m_tinker_mod(light, locy);
extern "C" int* m_tinker_mod(light, locz);
extern "C" int* m_tinker_mod(light, rgx);
extern "C" int* m_tinker_mod(light, rgy);
extern "C" int* m_tinker_mod(light, rgz);

int& nlight = m_tinker_mod(light, nlight);
int*& kbx = m_tinker_mod(light, kbx);
int*& kby = m_tinker_mod(light, kby);
int*& kbz = m_tinker_mod(light, kbz);
int*& kex = m_tinker_mod(light, kex);
int*& key = m_tinker_mod(light, key);
int*& kez = m_tinker_mod(light, kez);
int*& locx = m_tinker_mod(light, locx);
int*& locy = m_tinker_mod(light, locy);
int*& locz = m_tinker_mod(light, locz);
int*& rgx = m_tinker_mod(light, rgx);
int*& rgy = m_tinker_mod(light, rgy);
int*& rgz = m_tinker_mod(light, rgz);
#endif

} TINKER_NAMESPACE_END

#endif
