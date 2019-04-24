#ifndef TINKER_MOD_FILES_HH_
#define TINKER_MOD_FILES_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace files {
extern int& nprior;
extern int& ldir;
extern int& leng;
extern char (&filename)[240];
extern char (&outfile)[240];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(files, nprior);
extern "C" int m_tinker_mod(files, ldir);
extern "C" int m_tinker_mod(files, leng);
extern "C" char m_tinker_mod(files, filename)[240];
extern "C" char m_tinker_mod(files, outfile)[240];

int& nprior = m_tinker_mod(files, nprior);
int& ldir = m_tinker_mod(files, ldir);
int& leng = m_tinker_mod(files, leng);
char (&filename)[240] = m_tinker_mod(files, filename);
char (&outfile)[240] = m_tinker_mod(files, outfile);
#endif

} TINKER_NAMESPACE_END

#endif
