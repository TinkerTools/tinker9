#ifndef TINKER_MOD_SEQUEN_HH_
#define TINKER_MOD_SEQUEN_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace sequen {
using namespace sizes;

extern int& nseq;
extern int& nchain;
extern int (&ichain)[maxres][2];
extern int (&seqtyp)[maxres];
extern char (&chnnam)[maxres][1];
extern char (&seq)[maxres][3];
extern char (&chntyp)[maxres][7];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(sequen, nseq);
extern "C" int m_tinker_mod(sequen, nchain);
extern "C" int m_tinker_mod(sequen, ichain)[maxres][2];
extern "C" int m_tinker_mod(sequen, seqtyp)[maxres];
extern "C" char m_tinker_mod(sequen, chnnam)[maxres][1];
extern "C" char m_tinker_mod(sequen, seq)[maxres][3];
extern "C" char m_tinker_mod(sequen, chntyp)[maxres][7];

int& nseq = m_tinker_mod(sequen, nseq);
int& nchain = m_tinker_mod(sequen, nchain);
int (&ichain)[maxres][2] = m_tinker_mod(sequen, ichain);
int (&seqtyp)[maxres] = m_tinker_mod(sequen, seqtyp);
char (&chnnam)[maxres][1] = m_tinker_mod(sequen, chnnam);
char (&seq)[maxres][3] = m_tinker_mod(sequen, seq);
char (&chntyp)[maxres][7] = m_tinker_mod(sequen, chntyp);
#endif

} TINKER_NAMESPACE_END

#endif
