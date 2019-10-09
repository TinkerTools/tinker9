#pragma once
#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"
#include <tinker/detail/ktrtor.hh>


TINKER_NAMESPACE_BEGIN
// module bitor
TINKER_EXTERN int nbitor;
TINKER_EXTERN device_pointer<int, 5> ibitor;

// module tortor
TINKER_EXTERN int ntortor;
TINKER_EXTERN device_pointer<int, 3> itt;

// module ktrtor
TINKER_EXTERN device_pointer<int> tnx, tny; // of size maxntt
// of size (maxtgrd,maxntt) i.e. [maxntt][maxtgrd]
TINKER_EXTERN device_pointer<real, ktrtor::maxtgrd> ttx, tty;
// of size (maxtgrd2,maxntt) i.e. [maxntt][maxtgrd2]
TINKER_EXTERN device_pointer<real, ktrtor::maxtgrd2> tbf, tbx, tby, tbxy;

TINKER_EXTERN device_pointer<int> chkttor_ia_; // of size ntortor

TINKER_EXTERN real ttorunit;

TINKER_EXTERN energy_buffer ett;
TINKER_EXTERN virial_buffer vir_ett;

void etortor_data(rc_op op);

void etortor(int vers);
TINKER_NAMESPACE_END
