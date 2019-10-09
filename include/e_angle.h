#ifndef TINKER_E_ANGLE_H_
#define TINKER_E_ANGLE_H_

#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class eangle_t : int
{
   in_plane,
   harmonic,
   linear,
   fourier
};

// module angbnd
TINKER_EXTERN int nangle;
TINKER_EXTERN device_pointer<int, 4> iang;
TINKER_EXTERN device_pointer<real> ak, anat;

// module angpot
TINKER_EXTERN real angunit;
TINKER_EXTERN real cang, qang, pang, sang;
TINKER_EXTERN device_pointer<eangle_t> angtyp;

TINKER_EXTERN energy_buffer ea;
TINKER_EXTERN virial_buffer vir_ea;

void eangle_data(rc_op op);

void eangle(int vers);
TINKER_NAMESPACE_END

#endif
