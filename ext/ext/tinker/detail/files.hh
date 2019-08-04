#ifndef TINKER_MOD_FILES_HH_
#define TINKER_MOD_FILES_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace files {
extern int& nprior;
extern int& ldir;
extern int& leng;
extern char (&filename)[240];
extern char (&outfile)[240];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(files, nprior);
extern "C" int TINKER_MOD(files, ldir);
extern "C" int TINKER_MOD(files, leng);
extern "C" char TINKER_MOD(files, filename)[240];
extern "C" char TINKER_MOD(files, outfile)[240];

int& nprior = TINKER_MOD(files, nprior);
int& ldir = TINKER_MOD(files, ldir);
int& leng = TINKER_MOD(files, leng);
char (&filename)[240] = TINKER_MOD(files, filename);
char (&outfile)[240] = TINKER_MOD(files, outfile);
#endif
} TINKER_NAMESPACE_END

#endif
