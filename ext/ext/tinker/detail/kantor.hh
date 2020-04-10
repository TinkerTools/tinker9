#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kantor {
const int maxnat = 500;
extern double (&atcon)[maxnat][6];
extern char (&kat)[maxnat][16];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(kantor, atcon)[maxnat][6];
extern "C" char TINKER_MOD(kantor, kat)[maxnat][16];

double (&atcon)[maxnat][6] = TINKER_MOD(kantor, atcon);
char (&kat)[maxnat][16] = TINKER_MOD(kantor, kat);
#endif
} TINKER_NAMESPACE_END
