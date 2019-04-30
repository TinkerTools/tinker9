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
extern "C" int TINKER_MOD(inform, digits);
extern "C" int TINKER_MOD(inform, iprint);
extern "C" int TINKER_MOD(inform, iwrite);
extern "C" int TINKER_MOD(inform, isend);
extern "C" int TINKER_MOD(inform, silent);
extern "C" int TINKER_MOD(inform, verbose);
extern "C" int TINKER_MOD(inform, debug);
extern "C" int TINKER_MOD(inform, holdup);
extern "C" int TINKER_MOD(inform, abort);

int& digits = TINKER_MOD(inform, digits);
int& iprint = TINKER_MOD(inform, iprint);
int& iwrite = TINKER_MOD(inform, iwrite);
int& isend = TINKER_MOD(inform, isend);
int& silent = TINKER_MOD(inform, silent);
int& verbose = TINKER_MOD(inform, verbose);
int& debug = TINKER_MOD(inform, debug);
int& holdup = TINKER_MOD(inform, holdup);
int& abort = TINKER_MOD(inform, abort);
#endif
} TINKER_NAMESPACE_END

#endif
