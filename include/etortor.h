#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"
#include <tinker/detail/ktrtor.hh>


namespace tinker {
// module bitor
TINKER_EXTERN int nbitor;
TINKER_EXTERN pointer<int, 5> ibitor;

// module tortor
TINKER_EXTERN int ntortor;
TINKER_EXTERN pointer<int, 3> itt;

// module ktrtor
TINKER_EXTERN pointer<int> tnx, tny; // of size maxntt
// of size (maxtgrd,maxntt) i.e. [maxntt][maxtgrd]
TINKER_EXTERN pointer<real, ktrtor::maxtgrd> ttx, tty;
// of size (maxtgrd2,maxntt) i.e. [maxntt][maxtgrd2]
TINKER_EXTERN pointer<real, ktrtor::maxtgrd2> tbf, tbx, tby, tbxy;

TINKER_EXTERN pointer<int> chkttor_ia_; // of size ntortor

TINKER_EXTERN real ttorunit;

void etortor_data(rc_op op);

void etortor(int vers);
void etortor_acc(int);
}
