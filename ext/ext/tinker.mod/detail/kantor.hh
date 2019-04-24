#ifndef TINKER_MOD_KANTOR_HH_
#define TINKER_MOD_KANTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kantor {
const int maxnat = 500;
extern double (&atcon)[maxnat][6];
extern char (&kat)[maxnat][16];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kantor, atcon)[maxnat][6];
extern "C" char m_tinker_mod(kantor, kat)[maxnat][16];

double (&atcon)[maxnat][6] = m_tinker_mod(kantor, atcon);
char (&kat)[maxnat][16] = m_tinker_mod(kantor, kat);
#endif

} TINKER_NAMESPACE_END

#endif
