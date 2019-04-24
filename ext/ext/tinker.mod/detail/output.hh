#ifndef TINKER_MOD_OUTPUT_HH_
#define TINKER_MOD_OUTPUT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace output {
extern int& archive;
extern int& noversion;
extern int& overwrite;
extern char (&coordtype)[9];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(output, archive);
extern "C" int m_tinker_mod(output, noversion);
extern "C" int m_tinker_mod(output, overwrite);
extern "C" char m_tinker_mod(output, coordtype)[9];

int& archive = m_tinker_mod(output, archive);
int& noversion = m_tinker_mod(output, noversion);
int& overwrite = m_tinker_mod(output, overwrite);
char (&coordtype)[9] = m_tinker_mod(output, coordtype);
#endif

} TINKER_NAMESPACE_END

#endif
