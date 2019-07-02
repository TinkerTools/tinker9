#ifndef TINKER_GPU_E_ANGLE_H_
#define TINKER_GPU_E_ANGLE_H_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  angle_in_plane = 0x001,
  angle_harmonic = 0x002,
  angle_linear = 0x004,
  angle_fourier = 0x008
};

// module angbnd
extern int nangle;
extern int (*iang)[4];
extern real *ak, *anat;

// module angpot
extern real angunit;
extern real cang, qang, pang, sang;
extern int* angtyp;

extern real* ea;
extern real* vir_ea;

void eangle_data(rc_t rc);

void eangle(int vers);
}
TINKER_NAMESPACE_END

#endif
