#ifndef TINKER_E_TORTOR_H_
#define TINKER_E_TORTOR_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "ext/tinker/detail/ktrtor.hh"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
// module bitor
TINKER_EXTERN int nbitor;
TINKER_EXTERN DeviceVector<int, 5> ibitor_vec;

// module tortor
TINKER_EXTERN int ntortor;
TINKER_EXTERN DeviceVector<int, 3> itt_vec;

// module ktrtor
TINKER_EXTERN DeviceVector<int> tnx_vec, tny_vec; // of size maxntt
// of size (maxtgrd,maxntt) i.e. [maxntt][maxtgrd]
TINKER_EXTERN DeviceVector<real, ktrtor::maxtgrd> ttx_vec, tty_vec;
// of size (maxtgrd2,maxntt) i.e. [maxntt][maxtgrd2]
TINKER_EXTERN DeviceVector<real, ktrtor::maxtgrd2> tbf_vec, tbx_vec, tby_vec,
    tbxy_vec;

TINKER_EXTERN DeviceVector<int> chkttor_ia_vec_; // of size ntortor

TINKER_EXTERN real ttorunit;

TINKER_EXTERN BondedEnergy ett_handle;

void etortor_data(rc_op op);

void etortor(int vers);
TINKER_NAMESPACE_END

#endif
