#ifndef TINKER_MOD_RESDUE_HH_
#define TINKER_MOD_RESDUE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace resdue {
const int maxamino = 38;
const int maxnuc = 12;
extern int (&ntyp)[maxamino];
extern int (&catyp)[maxamino];
extern int (&ctyp)[maxamino];
extern int (&hntyp)[maxamino];
extern int (&otyp)[maxamino];
extern int (&hatyp)[maxamino];
extern int (&cbtyp)[maxamino];
extern int (&nntyp)[maxamino];
extern int (&cantyp)[maxamino];
extern int (&cntyp)[maxamino];
extern int (&hnntyp)[maxamino];
extern int (&ontyp)[maxamino];
extern int (&hantyp)[maxamino];
extern int (&nctyp)[maxamino];
extern int (&cactyp)[maxamino];
extern int (&cctyp)[maxamino];
extern int (&hnctyp)[maxamino];
extern int (&octyp)[maxamino];
extern int (&hactyp)[maxamino];
extern int (&o5typ)[maxnuc];
extern int (&c5typ)[maxnuc];
extern int (&h51typ)[maxnuc];
extern int (&h52typ)[maxnuc];
extern int (&c4typ)[maxnuc];
extern int (&h4typ)[maxnuc];
extern int (&o4typ)[maxnuc];
extern int (&c1typ)[maxnuc];
extern int (&h1typ)[maxnuc];
extern int (&c3typ)[maxnuc];
extern int (&h3typ)[maxnuc];
extern int (&c2typ)[maxnuc];
extern int (&h21typ)[maxnuc];
extern int (&o2typ)[maxnuc];
extern int (&h22typ)[maxnuc];
extern int (&o3typ)[maxnuc];
extern int (&ptyp)[maxnuc];
extern int (&optyp)[maxnuc];
extern int (&h5ttyp)[maxnuc];
extern int (&h3ttyp)[maxnuc];
extern char (&amino1)[maxamino][1];
extern char (&nuclz1)[maxnuc][1];
extern char (&amino)[maxamino][3];
extern char (&nuclz)[maxnuc][3];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(resdue, ntyp)[maxamino];
extern "C" int m_tinker_mod(resdue, catyp)[maxamino];
extern "C" int m_tinker_mod(resdue, ctyp)[maxamino];
extern "C" int m_tinker_mod(resdue, hntyp)[maxamino];
extern "C" int m_tinker_mod(resdue, otyp)[maxamino];
extern "C" int m_tinker_mod(resdue, hatyp)[maxamino];
extern "C" int m_tinker_mod(resdue, cbtyp)[maxamino];
extern "C" int m_tinker_mod(resdue, nntyp)[maxamino];
extern "C" int m_tinker_mod(resdue, cantyp)[maxamino];
extern "C" int m_tinker_mod(resdue, cntyp)[maxamino];
extern "C" int m_tinker_mod(resdue, hnntyp)[maxamino];
extern "C" int m_tinker_mod(resdue, ontyp)[maxamino];
extern "C" int m_tinker_mod(resdue, hantyp)[maxamino];
extern "C" int m_tinker_mod(resdue, nctyp)[maxamino];
extern "C" int m_tinker_mod(resdue, cactyp)[maxamino];
extern "C" int m_tinker_mod(resdue, cctyp)[maxamino];
extern "C" int m_tinker_mod(resdue, hnctyp)[maxamino];
extern "C" int m_tinker_mod(resdue, octyp)[maxamino];
extern "C" int m_tinker_mod(resdue, hactyp)[maxamino];
extern "C" int m_tinker_mod(resdue, o5typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, c5typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, h51typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, h52typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, c4typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, h4typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, o4typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, c1typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, h1typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, c3typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, h3typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, c2typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, h21typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, o2typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, h22typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, o3typ)[maxnuc];
extern "C" int m_tinker_mod(resdue, ptyp)[maxnuc];
extern "C" int m_tinker_mod(resdue, optyp)[maxnuc];
extern "C" int m_tinker_mod(resdue, h5ttyp)[maxnuc];
extern "C" int m_tinker_mod(resdue, h3ttyp)[maxnuc];
extern "C" char m_tinker_mod(resdue, amino1)[maxamino][1];
extern "C" char m_tinker_mod(resdue, nuclz1)[maxnuc][1];
extern "C" char m_tinker_mod(resdue, amino)[maxamino][3];
extern "C" char m_tinker_mod(resdue, nuclz)[maxnuc][3];

int (&ntyp)[maxamino] = m_tinker_mod(resdue, ntyp);
int (&catyp)[maxamino] = m_tinker_mod(resdue, catyp);
int (&ctyp)[maxamino] = m_tinker_mod(resdue, ctyp);
int (&hntyp)[maxamino] = m_tinker_mod(resdue, hntyp);
int (&otyp)[maxamino] = m_tinker_mod(resdue, otyp);
int (&hatyp)[maxamino] = m_tinker_mod(resdue, hatyp);
int (&cbtyp)[maxamino] = m_tinker_mod(resdue, cbtyp);
int (&nntyp)[maxamino] = m_tinker_mod(resdue, nntyp);
int (&cantyp)[maxamino] = m_tinker_mod(resdue, cantyp);
int (&cntyp)[maxamino] = m_tinker_mod(resdue, cntyp);
int (&hnntyp)[maxamino] = m_tinker_mod(resdue, hnntyp);
int (&ontyp)[maxamino] = m_tinker_mod(resdue, ontyp);
int (&hantyp)[maxamino] = m_tinker_mod(resdue, hantyp);
int (&nctyp)[maxamino] = m_tinker_mod(resdue, nctyp);
int (&cactyp)[maxamino] = m_tinker_mod(resdue, cactyp);
int (&cctyp)[maxamino] = m_tinker_mod(resdue, cctyp);
int (&hnctyp)[maxamino] = m_tinker_mod(resdue, hnctyp);
int (&octyp)[maxamino] = m_tinker_mod(resdue, octyp);
int (&hactyp)[maxamino] = m_tinker_mod(resdue, hactyp);
int (&o5typ)[maxnuc] = m_tinker_mod(resdue, o5typ);
int (&c5typ)[maxnuc] = m_tinker_mod(resdue, c5typ);
int (&h51typ)[maxnuc] = m_tinker_mod(resdue, h51typ);
int (&h52typ)[maxnuc] = m_tinker_mod(resdue, h52typ);
int (&c4typ)[maxnuc] = m_tinker_mod(resdue, c4typ);
int (&h4typ)[maxnuc] = m_tinker_mod(resdue, h4typ);
int (&o4typ)[maxnuc] = m_tinker_mod(resdue, o4typ);
int (&c1typ)[maxnuc] = m_tinker_mod(resdue, c1typ);
int (&h1typ)[maxnuc] = m_tinker_mod(resdue, h1typ);
int (&c3typ)[maxnuc] = m_tinker_mod(resdue, c3typ);
int (&h3typ)[maxnuc] = m_tinker_mod(resdue, h3typ);
int (&c2typ)[maxnuc] = m_tinker_mod(resdue, c2typ);
int (&h21typ)[maxnuc] = m_tinker_mod(resdue, h21typ);
int (&o2typ)[maxnuc] = m_tinker_mod(resdue, o2typ);
int (&h22typ)[maxnuc] = m_tinker_mod(resdue, h22typ);
int (&o3typ)[maxnuc] = m_tinker_mod(resdue, o3typ);
int (&ptyp)[maxnuc] = m_tinker_mod(resdue, ptyp);
int (&optyp)[maxnuc] = m_tinker_mod(resdue, optyp);
int (&h5ttyp)[maxnuc] = m_tinker_mod(resdue, h5ttyp);
int (&h3ttyp)[maxnuc] = m_tinker_mod(resdue, h3ttyp);
char (&amino1)[maxamino][1] = m_tinker_mod(resdue, amino1);
char (&nuclz1)[maxnuc][1] = m_tinker_mod(resdue, nuclz1);
char (&amino)[maxamino][3] = m_tinker_mod(resdue, amino);
char (&nuclz)[maxnuc][3] = m_tinker_mod(resdue, nuclz);
#endif

} TINKER_NAMESPACE_END

#endif
