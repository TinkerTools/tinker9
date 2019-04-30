#ifndef TINKER_MOD_POLGRP_HH_
#define TINKER_MOD_POLGRP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace polgrp {
const int maxp11 = 120;
const int maxp12 = 120;
const int maxp13 = 120;
const int maxp14 = 120;
extern int*& np11;
extern int*& np12;
extern int*& np13;
extern int*& np14;
extern int*& ip11;
extern int*& ip12;
extern int*& ip13;
extern int*& ip14;

#ifdef TINKER_MOD_CPP_
extern "C" int* TINKER_MOD(polgrp, np11);
extern "C" int* TINKER_MOD(polgrp, np12);
extern "C" int* TINKER_MOD(polgrp, np13);
extern "C" int* TINKER_MOD(polgrp, np14);
extern "C" int* TINKER_MOD(polgrp, ip11);
extern "C" int* TINKER_MOD(polgrp, ip12);
extern "C" int* TINKER_MOD(polgrp, ip13);
extern "C" int* TINKER_MOD(polgrp, ip14);

int*& np11 = TINKER_MOD(polgrp, np11);
int*& np12 = TINKER_MOD(polgrp, np12);
int*& np13 = TINKER_MOD(polgrp, np13);
int*& np14 = TINKER_MOD(polgrp, np14);
int*& ip11 = TINKER_MOD(polgrp, ip11);
int*& ip12 = TINKER_MOD(polgrp, ip12);
int*& ip13 = TINKER_MOD(polgrp, ip13);
int*& ip14 = TINKER_MOD(polgrp, ip14);
#endif
} TINKER_NAMESPACE_END

#endif
