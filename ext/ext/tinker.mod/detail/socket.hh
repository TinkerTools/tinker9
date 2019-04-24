#ifndef TINKER_MOD_SOCKET_HH_
#define TINKER_MOD_SOCKET_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace socket {
extern int& skttyp;
extern int& cstep;
extern double& cdt;
extern double& cenergy;
extern int& sktstart;
extern int& sktstop;
extern int& use_socket;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(socket, skttyp);
extern "C" int m_tinker_mod(socket, cstep);
extern "C" double m_tinker_mod(socket, cdt);
extern "C" double m_tinker_mod(socket, cenergy);
extern "C" int m_tinker_mod(socket, sktstart);
extern "C" int m_tinker_mod(socket, sktstop);
extern "C" int m_tinker_mod(socket, use_socket);

int& skttyp = m_tinker_mod(socket, skttyp);
int& cstep = m_tinker_mod(socket, cstep);
double& cdt = m_tinker_mod(socket, cdt);
double& cenergy = m_tinker_mod(socket, cenergy);
int& sktstart = m_tinker_mod(socket, sktstart);
int& sktstop = m_tinker_mod(socket, sktstop);
int& use_socket = m_tinker_mod(socket, use_socket);
#endif

} TINKER_NAMESPACE_END

#endif
