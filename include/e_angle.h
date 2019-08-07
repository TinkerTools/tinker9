#ifndef TINKER_E_ANGLE_H_
#define TINKER_E_ANGLE_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum {
  angle_in_plane = 0x001,
  angle_harmonic = 0x002,
  angle_linear = 0x004,
  angle_fourier = 0x008
};

// module angbnd
TINKER_EXTERN int nangle;
TINKER_EXTERN int (*iang)[4];
TINKER_EXTERN real *ak, *anat;

// module angpot
TINKER_EXTERN real angunit;
TINKER_EXTERN real cang, qang, pang, sang;
TINKER_EXTERN int* angtyp;

TINKER_EXTERN real* ea;
TINKER_EXTERN real* vir_ea;

void eangle_data(rc_op op);

void eangle(int vers);
TINKER_NAMESPACE_END

#endif
