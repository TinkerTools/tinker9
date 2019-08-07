#ifndef TINKER_E_TORTOR_H_
#define TINKER_E_TORTOR_H_

#include "rc_man.h"
#include <ext/tinker/detail/ktrtor.hh>

TINKER_NAMESPACE_BEGIN
// module bitor
TINKER_EXTERN int nbitor;
TINKER_EXTERN int (*ibitor)[5];

// module tortor
TINKER_EXTERN int ntortor;
TINKER_EXTERN int (*itt)[3];

// module ktrtor
TINKER_EXTERN int *tnx, *tny; // of size maxntt
// of size (maxtgrd,maxntt) i.e. [maxntt][maxtgrd]
TINKER_EXTERN real (*ttx)[ktrtor::maxtgrd];
TINKER_EXTERN real (*tty)[ktrtor::maxtgrd];
// of size (maxtgrd2,maxntt) i.e. [maxntt][maxtgrd2]
TINKER_EXTERN real (*tbf)[ktrtor::maxtgrd2];
TINKER_EXTERN real (*tbx)[ktrtor::maxtgrd2];
TINKER_EXTERN real (*tby)[ktrtor::maxtgrd2];
TINKER_EXTERN real (*tbxy)[ktrtor::maxtgrd2];

TINKER_EXTERN int* chkttor_ia_; // of size ntortor

TINKER_EXTERN real ttorunit;

TINKER_EXTERN real* ett;
TINKER_EXTERN real* vir_ett;

void etortor_data(rc_op op);

void etortor(int vers);
TINKER_NAMESPACE_END

#endif
