#ifndef TINKER_MOD_INFORM_HH_
#define TINKER_MOD_INFORM_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace inform {
const int maxask = 5;
extern int& digits;
extern int& iprint;
extern int& iwrite;
extern int& isend;
extern int& silent;
extern int& verbose;
extern int& debug;
extern int& holdup;
extern int& abort;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(inform, digits);
extern "C" int m_tinker_mod(inform, iprint);
extern "C" int m_tinker_mod(inform, iwrite);
extern "C" int m_tinker_mod(inform, isend);
extern "C" int m_tinker_mod(inform, silent);
extern "C" int m_tinker_mod(inform, verbose);
extern "C" int m_tinker_mod(inform, debug);
extern "C" int m_tinker_mod(inform, holdup);
extern "C" int m_tinker_mod(inform, abort);

int& digits = m_tinker_mod(inform, digits);
int& iprint = m_tinker_mod(inform, iprint);
int& iwrite = m_tinker_mod(inform, iwrite);
int& isend = m_tinker_mod(inform, isend);
int& silent = m_tinker_mod(inform, silent);
int& verbose = m_tinker_mod(inform, verbose);
int& debug = m_tinker_mod(inform, debug);
int& holdup = m_tinker_mod(inform, holdup);
int& abort = m_tinker_mod(inform, abort);
#endif

} TINKER_NAMESPACE_END

#endif
