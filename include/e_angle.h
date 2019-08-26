#ifndef TINKER_E_ANGLE_H_
#define TINKER_E_ANGLE_H_

#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class eangle_t : int { in_plane, harmonic, linear, fourier };

// module angbnd
TINKER_EXTERN int nangle;
TINKER_EXTERN int (*iang)[4];
TINKER_EXTERN real *ak, *anat;

// module angpot
TINKER_EXTERN real angunit;
TINKER_EXTERN real cang, qang, pang, sang;
TINKER_EXTERN eangle_t* angtyp;

TINKER_EXTERN BondedEnergy ea_handle;

void eangle_data(rc_op op);

void eangle(int vers);
TINKER_NAMESPACE_END

#endif
