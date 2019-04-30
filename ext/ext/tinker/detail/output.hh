#ifndef TINKER_MOD_OUTPUT_HH_
#define TINKER_MOD_OUTPUT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace output {
extern int& archive;
extern int& noversion;
extern int& overwrite;
extern char (&coordtype)[9];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(output, archive);
extern "C" int TINKER_MOD(output, noversion);
extern "C" int TINKER_MOD(output, overwrite);
extern "C" char TINKER_MOD(output, coordtype)[9];

int& archive = TINKER_MOD(output, archive);
int& noversion = TINKER_MOD(output, noversion);
int& overwrite = TINKER_MOD(output, overwrite);
char (&coordtype)[9] = TINKER_MOD(output, coordtype);
#endif
} TINKER_NAMESPACE_END

#endif
