#ifndef TINKER_MOD_CHUNKS_HH_
#define TINKER_MOD_CHUNKS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace chunks {
extern int& nchunk;
extern int& nchk1;
extern int& nchk2;
extern int& nchk3;
extern int& ngrd1;
extern int& ngrd2;
extern int& ngrd3;
extern int& nlpts;
extern int& nrpts;
extern int& grdoff;
extern int*& pmetable;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(chunks, nchunk);
extern "C" int m_tinker_mod(chunks, nchk1);
extern "C" int m_tinker_mod(chunks, nchk2);
extern "C" int m_tinker_mod(chunks, nchk3);
extern "C" int m_tinker_mod(chunks, ngrd1);
extern "C" int m_tinker_mod(chunks, ngrd2);
extern "C" int m_tinker_mod(chunks, ngrd3);
extern "C" int m_tinker_mod(chunks, nlpts);
extern "C" int m_tinker_mod(chunks, nrpts);
extern "C" int m_tinker_mod(chunks, grdoff);
extern "C" int* m_tinker_mod(chunks, pmetable);

int& nchunk = m_tinker_mod(chunks, nchunk);
int& nchk1 = m_tinker_mod(chunks, nchk1);
int& nchk2 = m_tinker_mod(chunks, nchk2);
int& nchk3 = m_tinker_mod(chunks, nchk3);
int& ngrd1 = m_tinker_mod(chunks, ngrd1);
int& ngrd2 = m_tinker_mod(chunks, ngrd2);
int& ngrd3 = m_tinker_mod(chunks, ngrd3);
int& nlpts = m_tinker_mod(chunks, nlpts);
int& nrpts = m_tinker_mod(chunks, nrpts);
int& grdoff = m_tinker_mod(chunks, grdoff);
int*& pmetable = m_tinker_mod(chunks, pmetable);
#endif

} TINKER_NAMESPACE_END

#endif
