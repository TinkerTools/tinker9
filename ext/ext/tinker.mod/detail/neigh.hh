#ifndef TINKER_MOD_NEIGH_HH_
#define TINKER_MOD_NEIGH_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace neigh {
extern int& maxvlst;
extern int& maxelst;
extern int& maxulst;
extern int*& nvlst;
extern int*& vlst;
extern int*& nelst;
extern int*& elst;
extern int*& nulst;
extern int*& ulst;
extern double& lbuffer;
extern double& pbuffer;
extern double& lbuf2;
extern double& pbuf2;
extern double& vbuf2;
extern double& vbufx;
extern double& dbuf2;
extern double& dbufx;
extern double& cbuf2;
extern double& cbufx;
extern double& mbuf2;
extern double& mbufx;
extern double& ubuf2;
extern double& ubufx;
extern double*& xvold;
extern double*& yvold;
extern double*& zvold;
extern double*& xeold;
extern double*& yeold;
extern double*& zeold;
extern double*& xuold;
extern double*& yuold;
extern double*& zuold;
extern int& dovlst;
extern int& dodlst;
extern int& doclst;
extern int& domlst;
extern int& doulst;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(neigh, maxvlst);
extern "C" int m_tinker_mod(neigh, maxelst);
extern "C" int m_tinker_mod(neigh, maxulst);
extern "C" int* m_tinker_mod(neigh, nvlst);
extern "C" int* m_tinker_mod(neigh, vlst);
extern "C" int* m_tinker_mod(neigh, nelst);
extern "C" int* m_tinker_mod(neigh, elst);
extern "C" int* m_tinker_mod(neigh, nulst);
extern "C" int* m_tinker_mod(neigh, ulst);
extern "C" double m_tinker_mod(neigh, lbuffer);
extern "C" double m_tinker_mod(neigh, pbuffer);
extern "C" double m_tinker_mod(neigh, lbuf2);
extern "C" double m_tinker_mod(neigh, pbuf2);
extern "C" double m_tinker_mod(neigh, vbuf2);
extern "C" double m_tinker_mod(neigh, vbufx);
extern "C" double m_tinker_mod(neigh, dbuf2);
extern "C" double m_tinker_mod(neigh, dbufx);
extern "C" double m_tinker_mod(neigh, cbuf2);
extern "C" double m_tinker_mod(neigh, cbufx);
extern "C" double m_tinker_mod(neigh, mbuf2);
extern "C" double m_tinker_mod(neigh, mbufx);
extern "C" double m_tinker_mod(neigh, ubuf2);
extern "C" double m_tinker_mod(neigh, ubufx);
extern "C" double* m_tinker_mod(neigh, xvold);
extern "C" double* m_tinker_mod(neigh, yvold);
extern "C" double* m_tinker_mod(neigh, zvold);
extern "C" double* m_tinker_mod(neigh, xeold);
extern "C" double* m_tinker_mod(neigh, yeold);
extern "C" double* m_tinker_mod(neigh, zeold);
extern "C" double* m_tinker_mod(neigh, xuold);
extern "C" double* m_tinker_mod(neigh, yuold);
extern "C" double* m_tinker_mod(neigh, zuold);
extern "C" int m_tinker_mod(neigh, dovlst);
extern "C" int m_tinker_mod(neigh, dodlst);
extern "C" int m_tinker_mod(neigh, doclst);
extern "C" int m_tinker_mod(neigh, domlst);
extern "C" int m_tinker_mod(neigh, doulst);

int& maxvlst = m_tinker_mod(neigh, maxvlst);
int& maxelst = m_tinker_mod(neigh, maxelst);
int& maxulst = m_tinker_mod(neigh, maxulst);
int*& nvlst = m_tinker_mod(neigh, nvlst);
int*& vlst = m_tinker_mod(neigh, vlst);
int*& nelst = m_tinker_mod(neigh, nelst);
int*& elst = m_tinker_mod(neigh, elst);
int*& nulst = m_tinker_mod(neigh, nulst);
int*& ulst = m_tinker_mod(neigh, ulst);
double& lbuffer = m_tinker_mod(neigh, lbuffer);
double& pbuffer = m_tinker_mod(neigh, pbuffer);
double& lbuf2 = m_tinker_mod(neigh, lbuf2);
double& pbuf2 = m_tinker_mod(neigh, pbuf2);
double& vbuf2 = m_tinker_mod(neigh, vbuf2);
double& vbufx = m_tinker_mod(neigh, vbufx);
double& dbuf2 = m_tinker_mod(neigh, dbuf2);
double& dbufx = m_tinker_mod(neigh, dbufx);
double& cbuf2 = m_tinker_mod(neigh, cbuf2);
double& cbufx = m_tinker_mod(neigh, cbufx);
double& mbuf2 = m_tinker_mod(neigh, mbuf2);
double& mbufx = m_tinker_mod(neigh, mbufx);
double& ubuf2 = m_tinker_mod(neigh, ubuf2);
double& ubufx = m_tinker_mod(neigh, ubufx);
double*& xvold = m_tinker_mod(neigh, xvold);
double*& yvold = m_tinker_mod(neigh, yvold);
double*& zvold = m_tinker_mod(neigh, zvold);
double*& xeold = m_tinker_mod(neigh, xeold);
double*& yeold = m_tinker_mod(neigh, yeold);
double*& zeold = m_tinker_mod(neigh, zeold);
double*& xuold = m_tinker_mod(neigh, xuold);
double*& yuold = m_tinker_mod(neigh, yuold);
double*& zuold = m_tinker_mod(neigh, zuold);
int& dovlst = m_tinker_mod(neigh, dovlst);
int& dodlst = m_tinker_mod(neigh, dodlst);
int& doclst = m_tinker_mod(neigh, doclst);
int& domlst = m_tinker_mod(neigh, domlst);
int& doulst = m_tinker_mod(neigh, doulst);
#endif

} TINKER_NAMESPACE_END

#endif
