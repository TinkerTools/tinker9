#ifndef TINKER_GPU_E_TORTOR_H_
#define TINKER_GPU_E_TORTOR_H_

#include "decl_basic.h"
#include <ext/tinker/detail/ktrtor.hh>

TINKER_NAMESPACE_BEGIN
namespace gpu {
// module bitor
extern int nbitor;
extern int (*ibitor)[5];

// module tortor
extern int ntortor;
extern int (*itt)[3];

// module ktrtor
extern int *tnx, *tny; // of size maxntt
// of size (maxtgrd,maxntt) i.e. [maxntt][maxtgrd]
extern real (*ttx)[ktrtor::maxtgrd];
extern real (*tty)[ktrtor::maxtgrd];
// of size (maxtgrd2,maxntt) i.e. [maxntt][maxtgrd2]
extern real (*tbf)[ktrtor::maxtgrd2];
extern real (*tbx)[ktrtor::maxtgrd2];
extern real (*tby)[ktrtor::maxtgrd2];
extern real (*tbxy)[ktrtor::maxtgrd2];

extern int* chkttor_ia__; // of size ntortor

extern real ttorunit;

extern real* ett;
extern real* vir_ett;

void etortor_data(rc_t rc);

void etortor_acc_impl__(int vers);
void etortor(int vers);
}
TINKER_NAMESPACE_END

#endif
