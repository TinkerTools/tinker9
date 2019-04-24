#ifndef TINKER_MOD_ACTION_HH_
#define TINKER_MOD_ACTION_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace action {
extern int& neb;
extern int& nea;
extern int& neba;
extern int& neub;
extern int& neaa;
extern int& neopb;
extern int& neopd;
extern int& neid;
extern int& neit;
extern int& net;
extern int& nept;
extern int& nebt;
extern int& neat;
extern int& nett;
extern int& nev;
extern int& ner;
extern int& nedsp;
extern int& nec;
extern int& necd;
extern int& ned;
extern int& nem;
extern int& nep;
extern int& nect;
extern int& new_;
extern int& nerxf;
extern int& nes;
extern int& nelf;
extern int& neg;
extern int& nex;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(action, neb);
extern "C" int m_tinker_mod(action, nea);
extern "C" int m_tinker_mod(action, neba);
extern "C" int m_tinker_mod(action, neub);
extern "C" int m_tinker_mod(action, neaa);
extern "C" int m_tinker_mod(action, neopb);
extern "C" int m_tinker_mod(action, neopd);
extern "C" int m_tinker_mod(action, neid);
extern "C" int m_tinker_mod(action, neit);
extern "C" int m_tinker_mod(action, net);
extern "C" int m_tinker_mod(action, nept);
extern "C" int m_tinker_mod(action, nebt);
extern "C" int m_tinker_mod(action, neat);
extern "C" int m_tinker_mod(action, nett);
extern "C" int m_tinker_mod(action, nev);
extern "C" int m_tinker_mod(action, ner);
extern "C" int m_tinker_mod(action, nedsp);
extern "C" int m_tinker_mod(action, nec);
extern "C" int m_tinker_mod(action, necd);
extern "C" int m_tinker_mod(action, ned);
extern "C" int m_tinker_mod(action, nem);
extern "C" int m_tinker_mod(action, nep);
extern "C" int m_tinker_mod(action, nect);
extern "C" int m_tinker_mod(action, new);
extern "C" int m_tinker_mod(action, nerxf);
extern "C" int m_tinker_mod(action, nes);
extern "C" int m_tinker_mod(action, nelf);
extern "C" int m_tinker_mod(action, neg);
extern "C" int m_tinker_mod(action, nex);

int& neb = m_tinker_mod(action, neb);
int& nea = m_tinker_mod(action, nea);
int& neba = m_tinker_mod(action, neba);
int& neub = m_tinker_mod(action, neub);
int& neaa = m_tinker_mod(action, neaa);
int& neopb = m_tinker_mod(action, neopb);
int& neopd = m_tinker_mod(action, neopd);
int& neid = m_tinker_mod(action, neid);
int& neit = m_tinker_mod(action, neit);
int& net = m_tinker_mod(action, net);
int& nept = m_tinker_mod(action, nept);
int& nebt = m_tinker_mod(action, nebt);
int& neat = m_tinker_mod(action, neat);
int& nett = m_tinker_mod(action, nett);
int& nev = m_tinker_mod(action, nev);
int& ner = m_tinker_mod(action, ner);
int& nedsp = m_tinker_mod(action, nedsp);
int& nec = m_tinker_mod(action, nec);
int& necd = m_tinker_mod(action, necd);
int& ned = m_tinker_mod(action, ned);
int& nem = m_tinker_mod(action, nem);
int& nep = m_tinker_mod(action, nep);
int& nect = m_tinker_mod(action, nect);
int& new_ = m_tinker_mod(action, new);
int& nerxf = m_tinker_mod(action, nerxf);
int& nes = m_tinker_mod(action, nes);
int& nelf = m_tinker_mod(action, nelf);
int& neg = m_tinker_mod(action, neg);
int& nex = m_tinker_mod(action, nex);
#endif

} TINKER_NAMESPACE_END

#endif
